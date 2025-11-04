#ifndef MD_LANGEVIN_EQUATION4_H
#define MD_LANGEVIN_EQUATION4_H

#include "array.h"
#include "STATISTICS.h"
#include "property.h"
#include "RANDOMNUMBER.h"
#include "simd512.h"
#include "simd.h"
#include "text.h"
#include "MOLCONST.h"
#include "streamer.h"
#include "clsHistogram.h"
#include "mathdef.h"
#include "DataLocator.h"
#include "cppheader.h"
#include "vector.h"
#include "clsTime.h"
#include "MOLUTILITY.h"
#include "inputParser.h"
#include "mpiKernel.h"

namespace LANGEVINEQUATION{
    void MeanSquareDisplacementBruteForce(std::vector<double> *r, std::vector<double>* result, std::vector<size_t>* intervalList, int first = -1, int last = -1){
        if( result->size() < intervalList->size() ){
            result->resize(intervalList->size());
        }
        if( first == -1 ) first = 0;
        if( first < 0 ) first = 0;
        if( last == -1 ) last = r->size();
        if( last > r->size() ) last = r->size();
        int size = last - first;

        for(int ii=0;ii<intervalList->size();ii++){
            int i = (*intervalList)[ii];
            (*result)[ii] = 0;
            if( last - i > 0 ){
                for(int k=first;k<last-i;k++){
                    double value = (*r)[k+i] - (*r)[k];
                    (*result)[ii] += value * value;
                }
                (*result)[ii] /= double(last-i-first);
            }
        }
    }
    void AutoCorrelationBruteForce(std::vector<double> *data, std::vector<double>* result, std::vector<size_t>* intervalList){
        for(size_t ii=0;ii<intervalList->size();ii++){
            size_t i = (*intervalList)[ii];
            (*result)[ii] = 0;
            for(size_t j=0;j<data->size()-i;j++){
                (*result)[ii] += (*data)[j] * (*data)[j + i];
            }
            (*result)[ii] = (*result)[ii] / (data->size() - i);
        }
    }

    class UserData{
    public:
        RANDOMNUMBER::UniformRandomNumberSequenceSIMD* m_RandomNumber;
        size_t m_iStep;
        bool m_bSameInitialBarrier;
        double m_bInitialValue[SIMDWIDTH];
        double m_EquilibriumConstant;
    };

    class BarrierStatus : public dataLocator::FiniteInifiniteIntegerTickArray1DElement{
    public:
        int m_State[SIMDWIDTH];
        size_t m_TimeStep[SIMDWIDTH];
    public:
        BarrierStatus() : FiniteInifiniteIntegerTickArray1DElement(){
            for(int i=0;i<SIMDWIDTH;i++){
                m_State[i] = 0;
                m_TimeStep[i] = 0;
            }
        }
    public:
        virtual void Initialize(int ID, void* userData){
            UserData* p = static_cast<UserData*>(userData);
            if( p->m_bSameInitialBarrier ){
                for(int i=0;i<SIMDWIDTH;i++){
                    m_State[i] = p->m_bInitialValue[i];
                    m_TimeStep[i] = p->m_iStep;
                }
            }else{
                __mSIMDd _prob = p->m_RandomNumber->GetRandomNumberSafe();
                alignas(SIMDBYTESIZE) double prob[SIMDWIDTH];
                _mmSIMD_store_pd(prob, _prob);

                double probState0 = 1.0 / (1.0 + p->m_EquilibriumConstant);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_State[i] = prob[i] < probState0 ? 0 : 1;
                    m_TimeStep[i] = p->m_iStep;
                }
            }
        };
        virtual std::string ToText(){
            std::string result = "";
            char buffer[256];
            for(int i=0;i<SIMDWIDTH;i++){
                sprintf(buffer, "%1d = %d, %10zd; ", i, m_State[i], m_TimeStep[i]);
                result += buffer;
            }
            return result;
        }
        virtual void Write(StreamerWriter *writer){
            for(int i=0;i<SIMDWIDTH;i++){
                writer->WriteInt(m_State[i]);
                writer->WriteSizeT(m_TimeStep[i]);
            }
        };
        virtual void Read(StreamerReader *reader){
            for(int i=0;i<SIMDWIDTH;i++){
                m_State[i] = reader->ReadInt();
                m_TimeStep[i] = reader->ReadSizeT();
            }
        };
    };

    class FirstPassageTimeCGData{
    public:
        std::vector<double> m_FirstPassageTimeCG2[SIMDWIDTH];
        double m_FirstPassageTimeLastTickCG2[SIMDWIDTH];
        bool m_FirstPassageTimeToUpdateCG2[SIMDWIDTH];
        int m_FirstPassageTimeInitialWellCG2[SIMDWIDTH];
    };

    class TrajectoryData{
    public:
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajX[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajV[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajF[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajPE[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajKE[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajEta[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajState[SIMDWIDTH];
        alignas(SIMDBYTESIZE) std::vector<double> m_TrajBarrier[SIMDWIDTH];
        __mSIMDd m_X;
        __mSIMDd m_V;
        __mSIMDd m_F;
        __mSIMDd m_E;
        __mSIMDd m_Eta;
        __mSIMDd m_Barrier;

        // for two state
        __mSIMDd m_CurrentState;
        __mSIMDd m_PreviousState;
        __mSIMDd m_SlopePositive;
        __mSIMDd m_SlopeNegative;
        __mSIMDd m_CurrentBarrierHeight;

        __mSIMDd m_StateCount[2];

        // for well and periodicity
        __mSIMDd m_WellIndex;
        __mSIMDd m_PeriodicIndex;
        __mSIMDd m_WellIndexPrevious;
        __mSIMDd m_PeriodicIndexPrevious;
        __mmask8 m_bOnEnteringDifferentWell;
        __mmask8 m_bOnEnteringDifferentPeriodicity;

        alignas(SIMDBYTESIZE) double m_V0Instance[SIMDWIDTH];
        alignas(SIMDBYTESIZE) double m_Eta0Instance[SIMDWIDTH];
        alignas(SIMDBYTESIZE) double m_BarrierState0Instance[SIMDWIDTH];

        alignas(SIMDBYTESIZE) double m_WellIndexStored[SIMDWIDTH];
        alignas(SIMDBYTESIZE) double m_PeriodicIndexStored[SIMDWIDTH];
        alignas(SIMDBYTESIZE) double m_WellIndexPreviousStored[SIMDWIDTH];
        alignas(SIMDBYTESIZE) double m_PeriodicIndexPreviousStored[SIMDWIDTH];
        bool m_bOnEnteringDifferentWellStored[SIMDWIDTH];
        bool m_bOnEnteringDifferentPeriodicityStored[SIMDWIDTH];

        // for waiting time 
        bool m_bWaitingTimeDistribution;
        alignas(SIMDBYTESIZE) double m_TimeTickPrevious[SIMDWIDTH];
        std::vector<double> m_BarrierTimeTick[SIMDWIDTH];
        bool m_BarrierPassed[SIMDWIDTH];

        // for first passage time
        bool m_bFirstPassageTime;
        // from bottom to top
        std::vector<double> m_FirstPassageTime[SIMDWIDTH];
        std::vector<double> m_InstantBarrierList[SIMDWIDTH];
        std::vector<double> m_InitialBarrierList[SIMDWIDTH];
        double m_FirstPassageTimeLastTick[SIMDWIDTH];
        bool m_FirstPassageTimeToUpdate[SIMDWIDTH];
        double m_InitialBarrierState[SIMDWIDTH][2];
        int m_InitialBarrierIndex[SIMDWIDTH][2];
        // from bottom to another bottom
        std::vector<double> m_FirstPassageTime2[SIMDWIDTH];
        double m_FirstPassageTimeLastTick2[SIMDWIDTH];
        bool m_FirstPassageTimeToUpdate2[SIMDWIDTH];
        int m_InitialBarrierIndex2[SIMDWIDTH][2];

        // for transitionPath time
        bool m_bTransitionPathTime;
        // from bottom to top
        std::vector<double> m_TransitionPathTime[SIMDWIDTH];
        double m_TransitionPathTimeLastTick[SIMDWIDTH];
        bool m_TransitionPathTimeToUpdate[SIMDWIDTH];
        // from bottom to another bottom
        std::vector<double> m_TransitionPathTime2[SIMDWIDTH];
        double m_TransitionPathTimeLastTick2[SIMDWIDTH];
        bool m_TransitionPathTimeToUpdate2[SIMDWIDTH];
        int m_TransitionPathTimeInitialWell2[SIMDWIDTH];

        // for first passage CG time
        bool m_bFirstPassageTimeCG;
        // from bottom to top, meaningless
        // std::vector<double> m_FirstPassageTimeCG[SIMDWIDTH];
        // double m_FirstPassageTimeLastTickCG[SIMDWIDTH];
        // bool m_FirstPassageTimeToUpdateCG[SIMDWIDTH];
        // from bottom to another bottom
        std::vector<int> m_FirstPassageTimeCGSize;
        std::vector<FirstPassageTimeCGData> m_FirstPassageTimeCGData;

        dataLocator::FiniteInifiniteIntegerTickArray1D<BarrierStatus> m_BarrierStatus;
    public:
        void Clear(){
            for(int i=0;i<SIMDWIDTH;i++){
                m_TrajX[i].clear();
                m_TrajV[i].clear();
                m_TrajF[i].clear();
                m_TrajPE[i].clear();
                m_TrajKE[i].clear();
                m_TrajEta[i].clear();
                m_TrajState[i].clear();
                m_TrajBarrier[i].clear();

                m_BarrierTimeTick[i].clear();
                m_FirstPassageTime[i].clear();
                m_InstantBarrierList[i].clear();
                m_InitialBarrierList[i].clear();

                m_FirstPassageTime2[i].clear();

                m_TransitionPathTime[i].clear();
                m_TransitionPathTime2[i].clear();

                // m_FirstPassageTimeCG[i].clear();
                for(int j=0;j<m_FirstPassageTimeCGData.size();j++){
                    m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i].clear();
                }
            }
            m_BarrierStatus.Release();
        }
        TrajectoryData(){
            for(int iChannel=0;iChannel<SIMDWIDTH;iChannel++){
                m_V0Instance[iChannel] = 0.0;
                m_Eta0Instance[iChannel] = 0.0;
                m_BarrierState0Instance[iChannel] = 0.0;

                m_WellIndexStored[iChannel] = 0.0;
                m_PeriodicIndexStored[iChannel] = 0.0;
                m_WellIndexPreviousStored[iChannel] = 0.0;
                m_PeriodicIndexPreviousStored[iChannel] = 0.0;
                m_bOnEnteringDifferentWellStored[iChannel] = false;
                m_bOnEnteringDifferentPeriodicityStored[iChannel] = false;
                
                m_TimeTickPrevious[iChannel] = 0.0;
                m_BarrierPassed[iChannel] = false;

                m_FirstPassageTimeLastTick[iChannel] = 0.0;
                m_FirstPassageTimeToUpdate[iChannel] = false;
                m_InitialBarrierState[iChannel][0] = 0.0;
                m_InitialBarrierState[iChannel][1] = 0.0;

                m_FirstPassageTimeLastTick2[iChannel] = 0.0;
                m_FirstPassageTimeToUpdate2[iChannel] = false;
                m_InitialBarrierIndex2[iChannel][0] = 0.0;
                m_InitialBarrierIndex2[iChannel][1] = 0.0;

                m_TransitionPathTimeLastTick[iChannel] = 0.0;
                m_TransitionPathTimeToUpdate[iChannel] = false;
                m_TransitionPathTimeLastTick2[iChannel] = 0.0;
                m_TransitionPathTimeToUpdate2[iChannel] = false;
                m_TransitionPathTimeInitialWell2[iChannel] = 0;
            }

            m_bWaitingTimeDistribution = false;
            m_bFirstPassageTime = false;

            m_bTransitionPathTime = false;

            m_X = _mmSIMD_set1_pd(0.0); 
            m_V = _mmSIMD_set1_pd(0.0); 
            m_E = _mmSIMD_set1_pd(0.0); 
            m_F = _mmSIMD_set1_pd(0.0); 
            m_Eta = _mmSIMD_set1_pd(0.0); 
            m_Barrier = _mmSIMD_set1_pd(0.0); 

            m_StateCount[0] = _mmSIMD_set1_pd(0.0);
            m_StateCount[1] = _mmSIMD_set1_pd(0.0);

            m_CurrentState = _mmSIMD_set1_pd(0.0); 
            m_PreviousState = _mmSIMD_set1_pd(0.0); 
            m_SlopePositive = _mmSIMD_set1_pd(0.0); 
            m_SlopeNegative = _mmSIMD_set1_pd(0.0); 
            m_CurrentBarrierHeight = _mmSIMD_set1_pd(0.0); 

            m_WellIndex = _mmSIMD_set1_pd(0.0); 
            m_PeriodicIndex = _mmSIMD_set1_pd(0.0); 
            m_WellIndexPrevious = _mmSIMD_set1_pd(0.0); 
            m_PeriodicIndexPrevious = _mmSIMD_set1_pd(0.0); 
            m_bOnEnteringDifferentWell = 0;
            m_bOnEnteringDifferentPeriodicity = 0;
        };
    public:
        virtual void Thin(size_t thinFreq, int ncore){
            for(int iChannel=0;iChannel<SIMDWIDTH;iChannel++){
                if( m_TrajX[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajX[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajX[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = m_TrajX[iChannel][f];
                    }
                    m_TrajX[iChannel] = buffer;
                }
                if( m_TrajV[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajV[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajV[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajV[iChannel][f]);
                    }
                    m_TrajV[iChannel] = buffer;
                }
                if( m_TrajF[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajF[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajF[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajF[iChannel][f]);
                    }
                    m_TrajF[iChannel] = buffer;
                }
                if( m_TrajPE[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajPE[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajPE[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajPE[iChannel][f]);
                    }
                    m_TrajPE[iChannel] = buffer;
                }
                if( m_TrajKE[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajKE[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajKE[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajKE[iChannel][f]);
                    }
                    m_TrajKE[iChannel] = buffer;
                }
                if( m_TrajEta[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajEta[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajEta[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajEta[iChannel][f]);
                    }
                    m_TrajEta[iChannel] = buffer;
                }
                if( m_TrajState[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajState[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajState[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajState[iChannel][f]);
                    }
                    m_TrajState[iChannel] = buffer;
                }
                if( m_TrajBarrier[iChannel].size() ){
                    std::vector<double> buffer;
                    buffer.resize((m_TrajBarrier[iChannel].size() - 1) / thinFreq + 1);
                    #pragma omp parallel for num_threads(ncore)
                    for(int f=0;f<m_TrajBarrier[iChannel].size();f+=thinFreq){
                        buffer[f/thinFreq] = (m_TrajBarrier[iChannel][f]);
                    }
                    m_TrajBarrier[iChannel] = buffer;
                }
            }
        };
        virtual void WriteRestartFile(StreamerWriter* writer){
            alignas(SIMDBYTESIZE) double x[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double v[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double f[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double energy[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double eta[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double barrier[SIMDWIDTH];

            alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double previousState[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double slopePositive[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double slopeNegative[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double currentBarrierHeight[SIMDWIDTH];

            alignas(SIMDBYTESIZE) double wellIndex[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double periodicIndex[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double wellIndexPrevious[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double periodicIndexPrevious[SIMDWIDTH];

            alignas(SIMDBYTESIZE) double stateCount[2][SIMDWIDTH];

            bool bOnEnteringDifferentWell[SIMDWIDTH];
            bool bOnEnteringDifferentPeriodicity[SIMDWIDTH];

            _mmSIMD_store_pd(x, m_X);
            _mmSIMD_store_pd(v, m_V);
            _mmSIMD_store_pd(f, m_F);
            _mmSIMD_store_pd(energy, m_E);
            _mmSIMD_store_pd(eta, m_Eta);
            _mmSIMD_store_pd(barrier, m_Barrier);

            _mmSIMD_store_pd(currentState, m_CurrentState);
            _mmSIMD_store_pd(previousState, m_PreviousState);
            _mmSIMD_store_pd(slopePositive, m_SlopePositive);
            _mmSIMD_store_pd(slopeNegative, m_SlopeNegative);
            _mmSIMD_store_pd(currentBarrierHeight, m_CurrentBarrierHeight);

            _mmSIMD_store_pd(wellIndex, m_WellIndex);
            _mmSIMD_store_pd(periodicIndex, m_PeriodicIndex);
            _mmSIMD_store_pd(wellIndexPrevious, m_WellIndexPrevious);
            _mmSIMD_store_pd(periodicIndexPrevious, m_PeriodicIndexPrevious);
            _mmask8_store_SIMD(bOnEnteringDifferentWell, m_bOnEnteringDifferentWell);
            _mmask8_store_SIMD(bOnEnteringDifferentPeriodicity, m_bOnEnteringDifferentPeriodicity);

            _mmSIMD_store_pd(stateCount[0], m_StateCount[0]);
            _mmSIMD_store_pd(stateCount[1], m_StateCount[1]);

            // write trajectory

            writer->WriteInt(9998);            
            writer->WriteInt(m_bWaitingTimeDistribution ? 1 : 0);
            writer->WriteInt(m_bFirstPassageTime ? 1 : 0);

            for(int i=0;i<SIMDWIDTH;i++){
                writer->WriteDouble(x[i]);
                writer->WriteDouble(v[i]);
                writer->WriteDouble(f[i]);
                writer->WriteDouble(energy[i]);
                writer->WriteDouble(eta[i]);
                writer->WriteDouble(barrier[i]);

                writer->WriteSizeT(m_TrajX[i].size());
                for(int j=0;j<m_TrajX[i].size();j++){
                    writer->WriteDouble(m_TrajX[i][j]);
                }

                writer->WriteSizeT(m_TrajV[i].size());
                for(int j=0;j<m_TrajV[i].size();j++){
                    writer->WriteDouble(m_TrajV[i][j]);
                }

                writer->WriteSizeT(m_TrajF[i].size());
                for(int j=0;j<m_TrajF[i].size();j++){
                    writer->WriteDouble(m_TrajF[i][j]);
                }

                writer->WriteSizeT(m_TrajPE[i].size());
                for(int j=0;j<m_TrajPE[i].size();j++){
                    writer->WriteDouble(m_TrajPE[i][j]);
                }

                writer->WriteSizeT(m_TrajKE[i].size());
                for(int j=0;j<m_TrajKE[i].size();j++){
                    writer->WriteDouble(m_TrajKE[i][j]);
                }

                writer->WriteSizeT(m_TrajEta[i].size());
                for(int j=0;j<m_TrajEta[i].size();j++){
                    writer->WriteDouble(m_TrajEta[i][j]);
                }

                writer->WriteSizeT(m_TrajState[i].size());
                for(int j=0;j<m_TrajState[i].size();j++){
                    writer->WriteDouble(m_TrajState[i][j]);
                }

                writer->WriteSizeT(m_TrajBarrier[i].size());
                for(int j=0;j<m_TrajBarrier[i].size();j++){
                    writer->WriteDouble(m_TrajBarrier[i][j]);
                }

                writer->WriteDouble(currentState[i]);
                writer->WriteDouble(previousState[i]);
                writer->WriteDouble(slopePositive[i]);
                writer->WriteDouble(slopeNegative[i]);
                writer->WriteDouble(currentBarrierHeight[i]);

                writer->WriteDouble(stateCount[0][i]);
                writer->WriteDouble(stateCount[1][i]);

                writer->WriteDouble(wellIndex[i]);
                writer->WriteDouble(periodicIndex[i]);
                writer->WriteDouble(wellIndexPrevious[i]);
                writer->WriteDouble(periodicIndexPrevious[i]);
                writer->WriteInt(bOnEnteringDifferentWell[i] ? 1 : 0);
                writer->WriteInt(bOnEnteringDifferentPeriodicity[i] ? 1 : 0);

                writer->WriteDouble(m_V0Instance[i]);
                writer->WriteDouble(m_Eta0Instance[i]);
                writer->WriteDouble(m_BarrierState0Instance[i]);

                writer->WriteDouble(m_WellIndexStored[i]);
                writer->WriteDouble(m_PeriodicIndexStored[i]);
                writer->WriteDouble(m_WellIndexPreviousStored[i]);
                writer->WriteDouble(m_PeriodicIndexPreviousStored[i]);
                writer->WriteInt(m_bOnEnteringDifferentWellStored[i] ? 1 : 0);
                writer->WriteInt(m_bOnEnteringDifferentPeriodicityStored[i] ? 1 : 0);

                if( m_bWaitingTimeDistribution ){
                    writer->WriteDouble(m_TimeTickPrevious[i]);
                    writer->WriteInt(m_BarrierPassed[i] ? 1 : 0);
                    writer->WriteSizeT(m_BarrierTimeTick[i].size());
                    writer->WriteVector(&m_BarrierTimeTick[i]);
                }

                if( m_bFirstPassageTime ){
                    writer->WriteSizeT(m_FirstPassageTime[i].size());
                    writer->WriteVector(&m_FirstPassageTime[i]);
                    writer->WriteSizeT(m_InstantBarrierList[i].size());
                    writer->WriteVector(&m_InstantBarrierList[i]);
                    writer->WriteSizeT(m_InitialBarrierList[i].size());
                    writer->WriteVector(&m_InitialBarrierList[i]);
                    writer->WriteDouble(m_FirstPassageTimeLastTick[i]);
                    writer->WriteInt(m_FirstPassageTimeToUpdate[i] ? 1 : 0);
                    writer->WriteDouble(m_InitialBarrierState[i][0]);
                    writer->WriteDouble(m_InitialBarrierState[i][1]);
                    writer->WriteInt(m_InitialBarrierIndex[i][0]);
                    writer->WriteInt(m_InitialBarrierIndex[i][1]);

                    writer->WriteSizeT(m_FirstPassageTime2[i].size());
                    writer->WriteVector(&m_FirstPassageTime2[i]);
                    writer->WriteDouble(m_FirstPassageTimeLastTick2[i]);
                    writer->WriteInt(m_FirstPassageTimeToUpdate2[i] ? 1 : 0);
                    writer->WriteInt(m_InitialBarrierIndex2[i][0]);
                    writer->WriteInt(m_InitialBarrierIndex2[i][1]);
                    
                }

                if( m_bTransitionPathTime ){
                    writer->WriteSizeT(m_TransitionPathTime[i].size());
                    writer->WriteVector(&m_TransitionPathTime[i]);
                    writer->WriteDouble(m_TransitionPathTimeLastTick[i]);
                    writer->WriteInt(m_TransitionPathTimeToUpdate[i] ? 1 : 0);

                    writer->WriteSizeT(m_TransitionPathTime2[i].size());
                    writer->WriteVector(&m_TransitionPathTime2[i]);
                    writer->WriteDouble(m_TransitionPathTimeLastTick2[i]);
                    writer->WriteInt(m_TransitionPathTimeToUpdate2[i] ? 1 : 0);
                    writer->WriteInt(m_TransitionPathTimeInitialWell2[i]);
                }

                if( m_bFirstPassageTimeCG ){
                    // writer->WriteSizeT(m_FirstPassageTimeCG[i].size());
                    // writer->WriteVector(&m_FirstPassageTimeCG[i]);
                    // writer->WriteDouble(m_FirstPassageTimeLastTickCG[i]);
                    // writer->WriteInt(m_FirstPassageTimeToUpdateCG[i] ? 1 : 0);
                    writer->WriteSizeT(m_FirstPassageTimeCGSize.size());
                    writer->WriteVector(&m_FirstPassageTimeCGSize);

                    for(int j=0;j<m_FirstPassageTimeCGData.size();j++){
                        writer->WriteSizeT(m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i].size());
                        writer->WriteVector(&m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i]);
                        writer->WriteDouble(m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i]);
                        writer->WriteInt(m_FirstPassageTimeCGData[j].m_FirstPassageTimeToUpdateCG2[i] ? 1 : 0);
                        writer->WriteInt(m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i]);
                    }
                }
            }

            m_BarrierStatus.Write(writer);

        }
        virtual void ReadRestartFile(StreamerReader* reader){
            Clear();

            alignas(SIMDBYTESIZE) double x[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double v[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double f[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double energy[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double eta[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double barrier[SIMDWIDTH];

            alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double previousState[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double slopePositive[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double slopeNegative[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double currentBarrierHeight[SIMDWIDTH];

            alignas(SIMDBYTESIZE) double wellIndex[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double periodicIndex[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double wellIndexPrevious[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double periodicIndexPrevious[SIMDWIDTH];

            alignas(SIMDBYTESIZE) double stateCount[2][SIMDWIDTH];

            bool bOnEnteringDifferentWell[SIMDWIDTH];
            bool bOnEnteringDifferentPeriodicity[SIMDWIDTH];

            // read trajectory
            if( reader->ReadInt() != 9998 ) fprintf(stderr, "Error: restart file corrupted\n");
            m_bWaitingTimeDistribution = reader->ReadInt() == 1;
            m_bFirstPassageTime = reader->ReadInt() == 1;

size_t sizeTraj = 0;
size_t sizeStatus = 0;
size_t sizeWaitingTime = 0;
size_t sizeFPT = 0;
size_t sizeTPT = 0;
size_t sizeFPTCG = 0;
size_t sizeBarrier = 0;
size_t sizePos;

            for(int i=0;i<SIMDWIDTH;i++){
                x[i] = reader->ReadDouble();
                v[i] = reader->ReadDouble();
                f[i] = reader->ReadDouble();
                energy[i] = reader->ReadDouble();
                eta[i] = reader->ReadDouble();
                barrier[i] = reader->ReadDouble();
sizePos = reader->GetPosition();

                size_t size;
                size = reader->ReadSizeT();
                m_TrajX[i].resize(size);
                for(int j=0;j<m_TrajX[i].size();j++){
                    m_TrajX[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajV[i].resize(size);
                for(int j=0;j<m_TrajV[i].size();j++){
                    m_TrajV[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajF[i].resize(size);
                for(int j=0;j<m_TrajF[i].size();j++){
                    m_TrajF[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajPE[i].resize(size);
                for(int j=0;j<m_TrajPE[i].size();j++){
                    m_TrajPE[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajKE[i].resize(size);
                for(int j=0;j<m_TrajKE[i].size();j++){
                    m_TrajKE[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajEta[i].resize(size);
                for(int j=0;j<m_TrajEta[i].size();j++){
                    m_TrajEta[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajState[i].resize(size);
                for(int j=0;j<m_TrajState[i].size();j++){
                    m_TrajState[i][j] = reader->ReadDouble();
                }

                size = reader->ReadSizeT();
                m_TrajBarrier[i].resize(size);
                for(int j=0;j<m_TrajBarrier[i].size();j++){
                    m_TrajBarrier[i][j] = reader->ReadDouble();
                }
sizeTraj += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

                currentState[i] = reader->ReadDouble();
                previousState[i] = reader->ReadDouble();
                slopePositive[i] = reader->ReadDouble();
                slopeNegative[i] = reader->ReadDouble();
                currentBarrierHeight[i] = reader->ReadDouble();

                stateCount[0][i] = reader->ReadDouble();
                stateCount[1][i] = reader->ReadDouble();

                wellIndex[i] = reader->ReadDouble();
                periodicIndex[i] = reader->ReadDouble();
                wellIndexPrevious[i] = reader->ReadDouble();
                periodicIndexPrevious[i] = reader->ReadDouble();
                bOnEnteringDifferentWell[i] = reader->ReadInt() == 1;
                bOnEnteringDifferentPeriodicity[i] = reader->ReadInt() == 1;

                m_V0Instance[i] = reader->ReadDouble();
                m_Eta0Instance[i] = reader->ReadDouble();
                m_BarrierState0Instance[i] = reader->ReadDouble();

                m_WellIndexStored[i] = reader->ReadDouble();
                m_PeriodicIndexStored[i] = reader->ReadDouble();
                m_WellIndexPreviousStored[i] = reader->ReadDouble();
                m_PeriodicIndexPreviousStored[i] = reader->ReadDouble();
                m_bOnEnteringDifferentWellStored[i] = reader->ReadInt() == 1;
                m_bOnEnteringDifferentPeriodicityStored[i] = reader->ReadInt() == 1;
sizeStatus += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

                if( m_bWaitingTimeDistribution ){
                    m_TimeTickPrevious[i] = reader->ReadDouble();
                    m_BarrierPassed[i] = reader->ReadInt() == 1;
                    size = reader->ReadSizeT();
                    m_BarrierTimeTick[i].resize(size);
                    reader->ReadVector(&m_BarrierTimeTick[i][0], size);
                }
sizeWaitingTime += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

                if( m_bFirstPassageTime ){
                    size = reader->ReadSizeT();
                    m_FirstPassageTime[i].resize(size);
                    reader->ReadVector<double>(&m_FirstPassageTime[i][0], size);
                    size = reader->ReadSizeT();
                    m_InstantBarrierList[i].resize(size);
                    reader->ReadVector<double>(&m_InstantBarrierList[i][0], size);
                    size = reader->ReadSizeT();
                    m_InitialBarrierList[i].resize(size);
                    reader->ReadVector<double>(&m_InitialBarrierList[i][0], size);
                    m_FirstPassageTimeLastTick[i] = reader->ReadDouble();
                    m_FirstPassageTimeToUpdate[i] = reader->ReadInt() == 1;
                    m_InitialBarrierState[i][0] = reader->ReadDouble();
                    m_InitialBarrierState[i][1] = reader->ReadDouble();
                    m_InitialBarrierIndex[i][0] = reader->ReadInt();
                    m_InitialBarrierIndex[i][1] = reader->ReadInt();

                    size = reader->ReadSizeT();
                    m_FirstPassageTime2[i].resize(size);
                    reader->ReadVector<double>(&m_FirstPassageTime2[i][0], size);
                    m_FirstPassageTimeLastTick2[i] = reader->ReadDouble();
                    m_FirstPassageTimeToUpdate2[i] = reader->ReadInt() == 1;
                    m_InitialBarrierIndex2[i][0] = reader->ReadInt();
                    m_InitialBarrierIndex2[i][1] = reader->ReadInt();
                }
sizeFPT += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

                if( m_bTransitionPathTime ){
                    size = reader->ReadSizeT();
                    m_TransitionPathTime[i].resize(size);
                    reader->ReadVector<double>(&m_TransitionPathTime[i][0], size);
                    m_TransitionPathTimeLastTick[i] = reader->ReadDouble();
                    m_TransitionPathTimeToUpdate[i] = reader->ReadInt() == 1;

                    size = reader->ReadSizeT();
                    m_TransitionPathTime2[i].resize(size);
                    reader->ReadVector<double>(&m_TransitionPathTime2[i][0], size);
                    m_TransitionPathTimeLastTick2[i] = reader->ReadDouble();
                    m_TransitionPathTimeToUpdate2[i] = reader->ReadInt() == 1;
                    m_TransitionPathTimeInitialWell2[i] = reader->ReadInt();
                }
sizeTPT += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

                if( m_bFirstPassageTimeCG ){
                    // size = reader->ReadSizeT();
                    // m_FirstPassageTimeCG[i].resize(size);
                    // reader->ReadVector<double>(&m_FirstPassageTimeCG[i][0], size);
                    // m_FirstPassageTimeLastTickCG[i] = reader->ReadDouble();
                    // m_FirstPassageTimeToUpdateCG[i] = reader->ReadInt() == 1;

                    size = reader->ReadSizeT();
                    m_FirstPassageTimeCGSize.resize(size);
                    reader->ReadVector<int>(&m_FirstPassageTimeCGSize[0], size);

                    m_FirstPassageTimeCGData.resize(size);
                    for(int j=0;j<m_FirstPassageTimeCGData.size();j++){
                        size = reader->ReadSizeT();
                        m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i].resize(size);
                        reader->ReadVector<double>(&m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i][0], size);
                        m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i] = reader->ReadDouble();
                        m_FirstPassageTimeCGData[j].m_FirstPassageTimeToUpdateCG2[i] = reader->ReadInt() == 1;
                        m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i] = reader->ReadInt();
                    }
                }
sizeFPTCG += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

            }
            m_BarrierStatus.Read(reader);
sizeBarrier += reader->GetPosition() - sizePos;
sizePos = reader->GetPosition();

// fprintf(stderr, "%20s %20s\n", "Traj", text::FromByte2Readable(sizeTraj).c_str());
// fprintf(stderr, "%20s %20s\n", "Status", text::FromByte2Readable(sizeStatus).c_str());
// fprintf(stderr, "%20s %20s\n", "WaitingTime", text::FromByte2Readable(sizeWaitingTime).c_str());
// fprintf(stderr, "%20s %20s\n", "FPT", text::FromByte2Readable(sizeFPT).c_str());
// fprintf(stderr, "%20s %20s\n", "TPT", text::FromByte2Readable(sizeTPT).c_str());
// fprintf(stderr, "%20s %20s\n", "FPTCG", text::FromByte2Readable(sizeFPTCG).c_str());
// fprintf(stderr, "%20s %20s\n", "Barrier", text::FromByte2Readable(sizeBarrier).c_str());

            m_X = _mmSIMD_load_pd(x);
            m_V = _mmSIMD_load_pd(v);
            m_F = _mmSIMD_load_pd(f);
            m_E = _mmSIMD_load_pd(energy);
            m_Eta = _mmSIMD_load_pd(eta);
            m_Barrier = _mmSIMD_load_pd(barrier);

            m_CurrentState = _mmSIMD_load_pd(currentState);
            m_PreviousState = _mmSIMD_load_pd(previousState);
            m_SlopePositive = _mmSIMD_load_pd(slopePositive);
            m_SlopeNegative = _mmSIMD_load_pd(slopeNegative);
            m_CurrentBarrierHeight = _mmSIMD_load_pd(currentBarrierHeight);

            m_WellIndex = _mmSIMD_load_pd(wellIndex);
            m_PeriodicIndex = _mmSIMD_load_pd(periodicIndex);
            m_WellIndexPrevious = _mmSIMD_load_pd(wellIndexPrevious);
            m_PeriodicIndexPrevious = _mmSIMD_load_pd(periodicIndexPrevious);
            m_bOnEnteringDifferentWell = _mmask8_load_SIMD(bOnEnteringDifferentWell);
            m_bOnEnteringDifferentPeriodicity = _mmask8_load_SIMD(bOnEnteringDifferentPeriodicity);

            m_StateCount[0] = _mmSIMD_load_pd(stateCount[0]);
            m_StateCount[1] = _mmSIMD_load_pd(stateCount[1]);
        }
    }; // TrajectoryData

    class LangevinEquationSIMDMultiple{
    public:
        static constexpr size_t MAXRANDOMNUMBERBANKSIZE = 40000; // 7.6 M each

        double m_X0;
        double m_V0;
        double m_Gamma0; 
        double m_kB0;
        double m_T0;
        double m_M0;
        double m_TimeStep0;
        size_t m_NumberOfStep;

        size_t m_StatisticsFreq;
        size_t m_SaveFreq;
        bool m_SaveX;
        bool m_SaveV;
        bool m_SaveF;
        bool m_SaveE;
        bool m_SaveEta;
        bool m_SaveState;
        bool m_SaveBarrier;

        bool m_RandomizeVelocity;
        bool m_RandomizeEta;
        bool m_RandomizeState0;

        // eta
        double m_PeriodicLength0;
        double m_BarrierHeight0[2]; 
        double m_InitialState0;
        double m_SlopeHalfLength0;
        double m_Q0; // for compatibility not used for now
        double m_Tau0; // k21 = 1 / m_Tau, rate from high to low, the larger rate between k21 and k12
        double m_Eta0; // for compatibility not used for now

        bool m_bWaitingTimeDistribution;
        bool m_bFirstPassageTime;
        bool m_bTransitionPathTime;
        bool m_bFirstPassageTimeCG;

        std::vector<int> m_FirstPassageTimeCGSize;

        bool m_ReadRestartFile;
        bool m_WriteRestartFile;
        text m_InputRestartFileName;
        text m_OutputRestartFileName;

        int m_Rank;

        bool m_bSameInitialBarrier;

        bool m_bEqualProbability;
        double m_EquilibriumConstant; // K = K12, 
        double m_ProbabilityOfFirstState0;
    protected:
        
        // size_t m_NumberOfRandomNumberStep0;
        // size_t m_NumberOfRandomNumberPerStep;
        // size_t m_NumberOfSaveGroup;

        // m_NumberOfRandomSeeding1 * m_NumberOfRandomSeedingSize1 = m_SaveFreq;
        // size_t m_NumberOfRandomSeeding1;
        // size_t m_NumberOfRandomSeedingSize1;
        // m_NumberOfRandomSeeding1 * m_NumberOfRandomSeedingSize1 = m_numberOfFrame / m_SaveFreq;
        // size_t m_NumberOfRandomSeeding2;
        // size_t m_NumberOfRandomSeedingSize2;

        TrajectoryData* m_TrajectoryData;
        // run time variables
        __mSIMDd m_kBInverse;
        __mSIMDd m_HalfTimeStep;
        __mSIMDd m_ExpMinusGammaTimeStep;
        __mSIMDd m_M;
        __mSIMDd m_MHalf;
        __mSIMDd m_MInverse;
        __mSIMDd m_T;
        __mSIMDd m_Gamma; 
        __mSIMDd m_kB;
        __mSIMDd m_TimeStep;
        __mSIMDd m_VarianceOfRandomForceIntegration;
        double m_TauInverse0;
        double m_MinusTwoTimeStepTauInverse0;
        __mSIMDd m_TauInverse;
        int m_RandomNumberSeed;
        
        // eta
        __mSIMDd m_TransitionProbability;
        __mSIMDd m_TransitionProbabilityAsymm[2];
        __mSIMDd m_StateBarrierHeight_0;
        __mSIMDd m_StateBarrierHeight_1;
        __mSIMDd m_StateBarrierSlope_0;
        __mSIMDd m_StateBarrierSlope_1;
        __mSIMDd m_SlopeHalfLength;
        __mSIMDd m_SlopeBoundLow;
        __mSIMDd m_SlopeBoundHigh;

        __mSIMDd m_PeriodicLength;
        // __mSIMDd m_BarrierHeight;
        // __mSIMDd m_Q;
        __mSIMDd m_Tau; // k21 = 1 / m_Tau, rate from high to low, the larger rate between k21 and k12
        __mSIMDd m_HalfPeriodicLength;
        __mSIMDd m_PeriodicLengthInverse;
        // __mSIMDd m_OrensteinUhlenbeckVariance;
        // __mSIMDd m_OrensteinUnlenbackExpMinusGammaTimeStep;
        
        double m_PeriodicLengthInverse0;
        double m_HalfPeriodicLength0;

        int m_Signature;

        size_t m_NumberOfStepDonePreviously;

        long long m_TimeOfGaussianRandomNumberSIMDRenewPreviously;
        size_t m_progressOfGaussianRandomNumberSIMDPreviously;

        long long m_TimeOfUniformRandomNumberSIMDRenewPreviously;
        size_t m_progressOfUniformRandomNumberSIMDPreviously;

        long long m_TimeOfUniformRandomNumberRenewPreviously;
        size_t m_progressOfUniformRandomNumberPreviously;

        UserData m_UserData;

        __mSIMDd m_ProbabilityOfFirstState;

#ifdef DEBUG
        std::vector<size_t> m_StepCounter;
#endif
    protected:
        RANDOMNUMBER::GaussianRandomNumberSequenceSIMD m_GaussianRandomNumberSequenceSIMD;
        RANDOMNUMBER::UniformRandomNumberSequenceSIMD m_UniformRandomNumberSequenceSIMD;
        RANDOMNUMBER::UniformRandomNumberSequence m_UniformRandomNumberSequence;

        const __mSIMDd m_MinusOne = _mmSIMD_set1_pd(-1.0);
        const __mSIMDd m_MinusTwo = _mmSIMD_set1_pd(-2.0);
        const __mSIMDd m_One = _mmSIMD_set1_pd(1.0);
        const __mSIMDd m_Zero = _mmSIMD_set1_pd(0.0);
        const __mSIMDd m_Half = _mmSIMD_set1_pd(0.5);
    public:
        LangevinEquationSIMDMultiple(){

            // m_NumberOfRandomNumberStep0 = 2;
            // m_NumberOfRandomNumberPerStep = 2;

            m_RandomNumberSeed = -1;

            m_X0 = 0.0;
            m_V0 = 1.0;
            m_Gamma0 = 1.0;
            m_kB0 = 1.0;
            m_T0 = 1.0;
            m_M0 = 1.0;

            m_TimeStep0 = 0.1;
            m_NumberOfStep = 10000;

            m_StatisticsFreq = 1;
            m_SaveFreq = 1;
            m_SaveX = false;
            m_SaveV = false;
            m_SaveF = false;
            m_SaveE = false;

            m_RandomizeVelocity = false;
            m_RandomizeEta = false;

            // eta
            m_PeriodicLength0 = 5.0;
            m_BarrierHeight0[0] = 2.0;
            m_BarrierHeight0[1] = 6.0;
            m_InitialState0 = 0;
            m_RandomizeState0 = true;
            m_SlopeHalfLength0 = m_PeriodicLength0 * 0.5;
            
            m_Tau0 = 20.0;

            m_SaveEta = false;
            m_SaveState = false;
            m_SaveBarrier = false;
            
            m_Rank = 0;

            m_ReadRestartFile = false;
            m_InputRestartFileName = "";
            m_WriteRestartFile = false;
            m_OutputRestartFileName = "";

            m_Signature = 2345;

            m_NumberOfStepDonePreviously = 0;
            m_TimeOfGaussianRandomNumberSIMDRenewPreviously = 0;
            m_progressOfGaussianRandomNumberSIMDPreviously = 0;

            m_TimeOfUniformRandomNumberSIMDRenewPreviously = 0;
            m_progressOfUniformRandomNumberSIMDPreviously = 0;

            m_TimeOfUniformRandomNumberRenewPreviously = 0;
            m_progressOfUniformRandomNumberPreviously = 0;

            m_bFirstPassageTime = false;
            m_bTransitionPathTime = false;
            m_bFirstPassageTimeCG = false;

            m_bSameInitialBarrier = false;

            m_bEqualProbability = true;
            m_EquilibriumConstant = 1.0;
            m_ProbabilityOfFirstState0 = 0.5;
        };
        ~LangevinEquationSIMDMultiple(){
            // SAFE_DELETE(m_Generator);
            // SAFE_DELETE(m_GaussianWhiteNoise);
            SAFE_DELETE(m_TrajectoryData);
        }
    protected:
        virtual void Allocate(){
            try{
                m_TrajectoryData = new TrajectoryData;
            }catch ( std::bad_alloc& e){
                fprintf(stderr, "Error: failed to allocate trajectory data %s\n", e.what());
                exit(0);
            }
        }
    public:
        virtual bool WriteRestartFile(int iTrial){
            StreamerWriter writer;
            char buffer[256];
            sprintf(buffer, "%s_%05d_%05d.restart", 
                m_OutputRestartFileName.c_str(), m_Rank, iTrial);
            if( writer.Open(buffer, std::ios::binary) ){
                writer.WriteInt(m_Signature);
                writer.WriteInt(SIMDWIDTH);
                // write settings
                writer.WriteDouble(m_X0);
                writer.WriteDouble(m_V0);
                writer.WriteDouble(m_Gamma0);
                writer.WriteDouble(m_kB0);
                writer.WriteDouble(m_T0);
                writer.WriteDouble(m_M0);

                writer.WriteDouble(m_TimeStep0);
                writer.WriteSizeT(m_NumberOfStepDonePreviously + m_NumberOfStep);

                writer.WriteInt(m_bWaitingTimeDistribution ? 1 : 0);
                writer.WriteSizeT(m_StatisticsFreq);
                writer.WriteSizeT(m_SaveFreq);
                writer.WriteInt(m_SaveX ? 1 : 0);
                writer.WriteInt(m_SaveV ? 1 : 0);
                writer.WriteInt(m_SaveF ? 1 : 0);
                writer.WriteInt(m_SaveE ? 1 : 0);

                writer.WriteInt(m_RandomizeVelocity ? 1 : 0);
                writer.WriteInt(m_RandomizeEta ? 1 : 0);

                writer.WriteDouble(m_PeriodicLength0);
                writer.WriteDouble(m_BarrierHeight0[0]);
                writer.WriteDouble(m_BarrierHeight0[1]);
                writer.WriteDouble(m_InitialState0);
                writer.WriteInt(m_RandomizeState0 ? 1 : 0);
                writer.WriteInt(m_bSameInitialBarrier ? 1 : 0);
                writer.WriteDouble(m_Q0);
                writer.WriteDouble(m_Tau0);
                writer.WriteDouble(m_Eta0);

                writer.WriteInt(m_SaveEta ? 1 : 0);
                writer.WriteInt(m_SaveState ? 1 : 0);
                writer.WriteInt(m_SaveBarrier ? 1 : 0);

                writer.WriteInt(m_Rank);

                writer.Write<long long>(m_GaussianRandomNumberSequenceSIMD.GetRenewTime());
                size_t progressOfRandomNumberSIMDPreviously = -1;
                if( progressOfRandomNumberSIMDPreviously == -1 ){
                    progressOfRandomNumberSIMDPreviously = m_GaussianRandomNumberSequenceSIMD.GetProcess();
                }else if( m_GaussianRandomNumberSequenceSIMD.GetProcess() != progressOfRandomNumberSIMDPreviously ){
                    fprintf(stderr, "Warning: inconsistent progressOfRandomNumberSIMDPreviously\n");
                }
                writer.WriteSizeT(progressOfRandomNumberSIMDPreviously);

                writer.Write<long long>(m_UniformRandomNumberSequenceSIMD.GetRenewTime());
                size_t progressOfUniformRandomNumberSIMDPreviously = -1;
                if( progressOfUniformRandomNumberSIMDPreviously == -1 ){
                    progressOfUniformRandomNumberSIMDPreviously = m_UniformRandomNumberSequenceSIMD.GetProcess();
                }else if( m_UniformRandomNumberSequenceSIMD.GetProcess() != progressOfUniformRandomNumberSIMDPreviously ){
                    fprintf(stderr, "Warning: inconsistent progressOfUniformRandomNumberSIMDPreviously\n");
                }
                writer.WriteSizeT(progressOfUniformRandomNumberSIMDPreviously);

                writer.Write<long long>(m_UniformRandomNumberSequence.GetRenewTime());
                size_t progressOfUniformRandomNumberPreviously = -1;
                if( progressOfUniformRandomNumberPreviously == -1 ){
                    progressOfUniformRandomNumberPreviously = m_UniformRandomNumberSequence.GetProcess();
                }else if( m_UniformRandomNumberSequence.GetProcess() != progressOfUniformRandomNumberPreviously ){
                    fprintf(stderr, "Warning: inconsistent progressOfUniformRandomNumberPreviously\n");
                }
                writer.WriteSizeT(progressOfUniformRandomNumberPreviously);

                // printf("%10lld %10zd\n\n", m_GaussianRandomNumberSequenceSIMD.GetRenewTime(), m_GaussianRandomNumberSequenceSIMD.GetProcess());
                // printf("%10lld %10zd\n\n", m_UniformRandomNumberSequenceSIMD.GetRenewTime(), m_UniformRandomNumberSequenceSIMD.GetProcess());
                // printf("%10lld %10zd\n\n", m_UniformRandomNumberSequence.GetRenewTime(), m_UniformRandomNumberSequence.GetProcess());

                // write trajectory data
                m_TrajectoryData->WriteRestartFile(&writer);

                writer.WriteInt(9990);
                writer.Close();
                return true;
            }else{
                fprintf(stderr, "Error: can't open restart file (%s)\n",
                    buffer);
                return false;
            }
        }
        virtual bool ReadRestartFile(int iTrial){
            StreamerReader reader;
            char buffer[256];
            sprintf(buffer, "%s_%05d_%05d.restart", 
                m_InputRestartFileName.c_str(), m_Rank, iTrial);
            if( reader.Open(buffer, std::ios::binary) ){
                if( reader.ReadInt() != m_Signature ){ fprintf(stderr, "Error: inconsistent signature\n"); reader.Close(); return false; }
                if( reader.ReadInt() != SIMDWIDTH ){ fprintf(stderr, "Error: inconsistent simd witdth\n"); reader.Close(); return false; }
                // write settings
                if( reader.ReadDouble() != m_X0 ){ fprintf(stderr, "Error: inconsistent x0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_V0 ){ fprintf(stderr, "Error: inconsistent v0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_Gamma0 ){ fprintf(stderr, "Error: inconsistent gamma0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_kB0){ fprintf(stderr, "Error: inconsistent kb0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_T0){ fprintf(stderr, "Error: inconsistent T0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_M0){ fprintf(stderr, "Error: inconsistent m0\n"); reader.Close(); return false; }

                if( reader.ReadDouble() != m_TimeStep0 ){ fprintf(stderr, "Error: inconsistent time step0\n"); reader.Close(); return false; }
                m_NumberOfStepDonePreviously = reader.ReadSizeT();

                if( reader.ReadInt() != (m_bWaitingTimeDistribution ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent bWaitingTimeDistribution\n"); reader.Close(); return false;}
                if( reader.ReadSizeT() != m_StatisticsFreq ){ fprintf(stderr, "Error: inconsistent statisticsFreq\n"); reader.Close(); return false; };
                if( reader.ReadSizeT() != m_SaveFreq ){ fprintf(stderr, "Error: inconsistent save freq\n"); reader.Close(); return false; }
                if( reader.ReadInt() != (m_SaveX ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent bSaveX\n"); reader.Close(); return false; }
                reader.ReadInt(); // bSaveV not to check this because it is modified afterwards
                if( reader.ReadInt() != (m_SaveF ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent bSaveF\n"); reader.Close(); return false; }
                if( reader.ReadInt() != (m_SaveE ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent bSaveE\n"); reader.Close(); return false; }

                if( reader.ReadInt() != (m_RandomizeVelocity ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent randomize velocity\n"); reader.Close(); return false; };
                if( reader.ReadInt() != (m_RandomizeEta ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent randomize eta\n"); reader.Close(); return false; };

                if( reader.ReadDouble() != m_PeriodicLength0 ){ fprintf(stderr, "Error: inconsistent periodic length0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_BarrierHeight0[0] ){ fprintf(stderr, "Error: inconsistent barrier height0 [0]\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_BarrierHeight0[1] ){ fprintf(stderr, "Error: inconsistent barrier height0 [1]\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_InitialState0 ){ fprintf(stderr, "Error: inconsistent initialState0\n"); reader.Close(); return false; }
                if( reader.ReadInt() != m_RandomizeState0 ? 1 : 0 ){ fprintf(stderr, "Error: inconsistent randomizeState0\n"); reader.Close(); return false; }
                if( reader.ReadInt() != m_bSameInitialBarrier ? 1 : 0 ){ fprintf(stderr, "Error: inconsistent bSameInitialBarrier\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_Q0 ){ fprintf(stderr, "Error: inconsistent q0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_Tau0 ){ fprintf(stderr, "Error: inconsistent tau0\n"); reader.Close(); return false; }
                if( reader.ReadDouble() != m_Eta0 ){ fprintf(stderr, "Error: inconsistent eta0\n"); reader.Close(); return false; }

                if( reader.ReadInt() != (m_SaveEta ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent bSaveEta\n"); reader.Close(); return false; }
                if( reader.ReadInt() != (m_SaveState ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent bSaveState\n"); reader.Close(); return false; }
                if( reader.ReadInt() != (m_SaveBarrier ? 1 : 0) ){ fprintf(stderr, "Error: inconsistent save barrier\n"); reader.Close(); return false; }

                if( reader.ReadInt() != m_Rank ){ fprintf(stderr, "Error: inconsistent rank\n"); reader.Close(); return false;}

                m_TimeOfGaussianRandomNumberSIMDRenewPreviously = reader.Read<long long>();
                m_progressOfGaussianRandomNumberSIMDPreviously = reader.ReadSizeT();

                m_TimeOfUniformRandomNumberSIMDRenewPreviously = reader.Read<long long>();
                m_progressOfUniformRandomNumberSIMDPreviously = reader.ReadSizeT();

                m_TimeOfUniformRandomNumberRenewPreviously = reader.Read<long long>();
                m_progressOfUniformRandomNumberPreviously = reader.ReadSizeT();

                // read trajectory data
                m_TrajectoryData->ReadRestartFile(&reader);

                if( reader.ReadInt() != 9990 ) fprintf(stderr, "Error: restart file corrupted\n");
                reader.Close();
                return true;
            }else{
                fprintf(stderr, "Error: can't open restart file (%s)\n",
                    buffer);
                return false;
            }
        }
        virtual bool ReadRestartFileAll(int iTrial, bool bRun = true){
            StreamerReader reader;
            char buffer[256];
            sprintf(buffer, "%s_%05d_%05d.restart", 
                m_InputRestartFileName.c_str(), m_Rank, iTrial);
            if( reader.Open(buffer, std::ios::binary) ){
                m_Signature = reader.ReadInt();
                if( reader.ReadInt() != SIMDWIDTH ){ fprintf(stderr, "Error: inconsistent simd witdth\n"); reader.Close(); return false; }
                // write settings
                m_X0 = reader.ReadDouble();
                m_V0 = reader.ReadDouble();
                m_Gamma0 = reader.ReadDouble();
                m_kB0 = reader.ReadDouble();
                m_T0 = reader.ReadDouble();
                m_M0 = reader.ReadDouble();

                m_TimeStep0 = reader.ReadDouble();
                m_NumberOfStep = reader.ReadSizeT();
                m_NumberOfStepDonePreviously = 0;

                m_bWaitingTimeDistribution = reader.ReadInt() == 1;
                m_StatisticsFreq = reader.ReadSizeT();
                m_SaveFreq = reader.ReadSizeT();
                m_SaveX = reader.ReadInt() == 1;
                m_SaveV = reader.ReadInt() == 1;
                m_SaveF = reader.ReadInt() == 1;
                m_SaveE = reader.ReadInt() == 1;

                m_RandomizeVelocity = reader.ReadInt() == 1;
                m_RandomizeEta = reader.ReadInt() == 1;

                m_PeriodicLength0 = reader.ReadDouble();
                m_BarrierHeight0[0] = reader.ReadDouble();
                m_BarrierHeight0[1] = reader.ReadDouble();
                m_InitialState0 = reader.ReadDouble();
                m_RandomizeState0 = reader.ReadInt() == 1;
                m_bSameInitialBarrier = reader.ReadInt() == 1;
                m_Q0 = reader.ReadDouble();
                m_Tau0 = reader.ReadDouble();
                m_Eta0 = reader.ReadDouble();

                m_SaveEta = reader.ReadInt() == 1;
                m_SaveState = reader.ReadInt() == 1;
                m_SaveBarrier = reader.ReadInt() == 1;

                m_Rank = reader.ReadInt();

                m_TimeOfGaussianRandomNumberSIMDRenewPreviously = reader.Read<long long>();
                m_progressOfGaussianRandomNumberSIMDPreviously = reader.ReadSizeT();
                if( bRun ) m_GaussianRandomNumberSequenceSIMD.SetRenewTime(m_TimeOfGaussianRandomNumberSIMDRenewPreviously);

                m_TimeOfUniformRandomNumberSIMDRenewPreviously = reader.Read<long long>();
                m_progressOfUniformRandomNumberSIMDPreviously = reader.ReadSizeT();
                if( bRun ) m_UniformRandomNumberSequenceSIMD.SetRenewTime(m_TimeOfUniformRandomNumberSIMDRenewPreviously);

                m_TimeOfUniformRandomNumberRenewPreviously = reader.Read<long long>();
                m_progressOfUniformRandomNumberPreviously = reader.ReadSizeT();
                if( bRun ) m_UniformRandomNumberSequence.SetRenewTime(m_TimeOfUniformRandomNumberRenewPreviously);

                // read trajectory data
                m_TrajectoryData->ReadRestartFile(&reader);

                if( reader.ReadInt() != 9990 ) fprintf(stderr, "Error: restart file corrupted\n");

                reader.Close();
                return true;
            }else{
                fprintf(stderr, "Error: can't open restart file (%s)\n",
                    buffer);
                return false;
            }
        }
        virtual void Initialize(int randomNumberSeed){

            Allocate();
            m_TrajectoryData->m_bFirstPassageTime = m_bFirstPassageTime;
            m_TrajectoryData->m_bFirstPassageTimeCG = m_bFirstPassageTimeCG;
            m_TrajectoryData->m_bTransitionPathTime = m_bTransitionPathTime;
            m_TrajectoryData->m_bWaitingTimeDistribution = m_bWaitingTimeDistribution;
            m_TrajectoryData->m_FirstPassageTimeCGSize = m_FirstPassageTimeCGSize;
            m_TrajectoryData->m_FirstPassageTimeCGData.resize(m_FirstPassageTimeCGSize.size());

            m_RandomNumberSeed = randomNumberSeed;

            m_GaussianRandomNumberSequenceSIMD.Initialize(MAXRANDOMNUMBERBANKSIZE, m_RandomNumberSeed, 0.0, 1.0);
            m_UniformRandomNumberSequenceSIMD.Initialize(MAXRANDOMNUMBERBANKSIZE, m_RandomNumberSeed, 0.0, 1.0);
            m_UniformRandomNumberSequence.Initialize(MAXRANDOMNUMBERBANKSIZE, m_RandomNumberSeed, 0.0, 1.0);
        }
        virtual void EnergyAndForce(){
            // default is (x^2-1)^4, x belongs to [-1, 1]

            __mSIMDd _reducedX = _mmSIMD_mul_pd(m_TrajectoryData->m_X, m_PeriodicLengthInverse);

            __mSIMDd _periodicIndex = _mmSIMD_round_pd(_reducedX, _MM_FROUND_TO_NEAREST_INT);

            __mSIMDd _xp = _mmSIMD_sub_pd(
                m_TrajectoryData->m_X, _mmSIMD_mul_pd(_periodicIndex, m_PeriodicLength));

            __mmask8 _bOnNegative = _mmask8_and_SIMD(
                _mmSIMD_cmp_pd_mask(_xp, m_Zero, _CMP_LT_OS),            // x < 0
                _mmSIMD_cmp_pd_mask(_xp, m_SlopeBoundLow, _CMP_GT_OS)); // x > low

            __mmask8 _bOnPositive = _mmask8_and_SIMD(
                _mmSIMD_cmp_pd_mask(_xp, m_SlopeBoundHigh, _CMP_LT_OS),            // x < high
                _mmSIMD_cmp_pd_mask(_xp, m_Zero, _CMP_GT_OS)); // x > 0

            __mSIMDd _energy = m_Zero;
            //                              sub because slope is negative slope
            __mSIMDd _valueNegative = _mmSIMD_sub_pd(m_TrajectoryData->m_CurrentBarrierHeight, _mmSIMD_mul_pd(_xp, m_TrajectoryData->m_SlopeNegative));
            __mSIMDd _valuePositive = _mmSIMD_sub_pd(m_TrajectoryData->m_CurrentBarrierHeight, _mmSIMD_mul_pd(_xp, m_TrajectoryData->m_SlopePositive));
            _energy = _mmSIMD_mask_blend_pd(_bOnNegative, m_Zero, _valueNegative);
            _energy = _mmSIMD_mask_blend_pd(_bOnPositive, _energy, _valuePositive);

            __mSIMDd _force = m_Zero;
            _force = _mmSIMD_mask_blend_pd(_bOnNegative, m_Zero, m_TrajectoryData->m_SlopeNegative);
            _force = _mmSIMD_mask_blend_pd(_bOnPositive, _force, m_TrajectoryData->m_SlopePositive);

            m_TrajectoryData->m_E = _energy;
            m_TrajectoryData->m_F = _force;
            m_TrajectoryData->m_Barrier = m_TrajectoryData->m_CurrentBarrierHeight;
        };
        __mSIMDd KineticEnergy(__mSIMDd v){
            __mSIMDd _value = _mmSIMD_mul_pd(v, v);
            _value = _mmSIMD_mul_pd(m_MHalf, _value);
            return _value;
        }
        virtual bool PreIntegrate(int iTrial){
            // restart variables
            // m_GaussianRandomNumberSequenceSIMD.ResetRenewTime();
            m_TimeOfGaussianRandomNumberSIMDRenewPreviously = m_GaussianRandomNumberSequenceSIMD.GetRenewTime();
            m_progressOfGaussianRandomNumberSIMDPreviously = m_GaussianRandomNumberSequenceSIMD.GetProcess();

            // m_UniformRandomNumberSequenceSIMD.ResetRenewTime();
            m_TimeOfUniformRandomNumberSIMDRenewPreviously = m_UniformRandomNumberSequenceSIMD.GetRenewTime();
            m_progressOfUniformRandomNumberSIMDPreviously = m_UniformRandomNumberSequenceSIMD.GetProcess();

            // m_UniformRandomNumberSequence.ResetRenewTime();
            m_TimeOfUniformRandomNumberRenewPreviously = m_UniformRandomNumberSequence.GetRenewTime();
            m_progressOfUniformRandomNumberPreviously = m_UniformRandomNumberSequence.GetProcess();

            if( m_ReadRestartFile ){
                if( !ReadRestartFile(iTrial) ){
                    fprintf(stderr, "Error: can't read restart file on rank %10d\n", m_Rank);
                    return false;
                }
                // m_GaussianRandomNumberSequenceSIMD.SetRenewTime(m_TimeOfGaussianRandomNumberSIMDRenewPreviously);
                // m_UniformRandomNumberSequenceSIMD.SetRenewTime(m_TimeOfUniformRandomNumberSIMDRenewPreviously);
                // m_UniformRandomNumberSequence.SetRenewTime(m_TimeOfUniformRandomNumberRenewPreviously);

                // m_GaussianRandomNumberSequenceSIMD.SetProcess(m_progressOfGaussianRandomNumberSIMDPreviously);
                // m_UniformRandomNumberSequenceSIMD.SetProcess(m_progressOfUniformRandomNumberSIMDPreviously);
                // m_UniformRandomNumberSequence.SetProcess(m_progressOfUniformRandomNumberPreviously);

            // skip random numebr if necessary
            // if( m_TimeOfGaussianRandomNumberSIMDRenewPreviously != 0 ||
            //     m_progressOfGaussianRandomNumberSIMDPreviously != 0 ){
                m_GaussianRandomNumberSequenceSIMD.AdvanceFromZero(
                    m_TimeOfGaussianRandomNumberSIMDRenewPreviously, m_progressOfGaussianRandomNumberSIMDPreviously);
            // }
            // if( m_TimeOfUniformRandomNumberSIMDRenewPreviously != 0 ||
            //     m_progressOfUniformRandomNumberSIMDPreviously != 0 ){
                m_UniformRandomNumberSequenceSIMD.AdvanceFromZero(
                    m_TimeOfUniformRandomNumberSIMDRenewPreviously, m_progressOfUniformRandomNumberSIMDPreviously);
            // }
            // if( m_TimeOfUniformRandomNumberRenewPreviously != 0 ||
            //     m_progressOfUniformRandomNumberPreviously != 0 ){
                m_UniformRandomNumberSequence.AdvanceFromZero(
                    m_TimeOfUniformRandomNumberRenewPreviously, m_progressOfUniformRandomNumberPreviously);
            // }
            }
            // initialize barrier
            m_UserData.m_iStep = 0;
            m_UserData.m_RandomNumber = &m_UniformRandomNumberSequenceSIMD;
            m_UserData.m_bSameInitialBarrier = m_bSameInitialBarrier;
            for(int i=0;i<SIMDWIDTH;i++){
                m_UserData.m_bInitialValue[i] = m_TrajectoryData->m_BarrierState0Instance[i];
            }
            m_UserData.m_EquilibriumConstant = m_EquilibriumConstant;
            m_TrajectoryData->m_BarrierStatus.SetUserData(&m_UserData);

            return true;
        }
        virtual bool Integrate(int iTrial, bool bTimer = false){
            PrepareRecord();

            if( !PreIntegrate(iTrial) ){
                return false;
            }

            IntegrateInitializeVariables();

            {
                size_t iStep = 0;
                if( m_NumberOfStepDonePreviously == 0 ){
                    IntegrateStep0();
                    Record(iStep);
                }
            }

            if( m_SaveFreq >= 10 ){
                LoopTimer loopTimer;
                if( bTimer ){
                    loopTimer.Start();
                    fprintf(stderr, "\n");
                }
                size_t ngroup = m_NumberOfStep / m_SaveFreq;
                if( m_Tau0 == 0.0 ){
                    for(size_t igroup=0;igroup<ngroup;igroup++){
                        size_t iStepBase = igroup * m_SaveFreq + 1;

                        if( bTimer ){
                            if( igroup == 0 || loopTimer.Tick(1.0) ){
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStepBase, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }
                        for(size_t ii=0;ii<m_SaveFreq;ii++){
                            size_t iStep = m_NumberOfStepDonePreviously + iStepBase + ii;
                            IntegrateStep(iStep);
                            /*if( iStep % m_StatisticsFreq == 0 )*/ UpdateStatistics(iStep);
                        }
                        Record(m_NumberOfStepDonePreviously + iStepBase + m_SaveFreq - 1);

                        if( bTimer ){
                            if( igroup == ngroup - 1 || loopTimer.Tick(1.0) ){
                                loopTimer.Update(iStepBase + m_SaveFreq - 1, m_NumberOfStep + 1);
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStepBase + m_SaveFreq, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }
                    }
                }else{
                    for(size_t igroup=0;igroup<ngroup;igroup++){
                        size_t iStepBase = igroup * m_SaveFreq + 1;

                        if( bTimer ){
                            if( igroup == 0 || loopTimer.Tick(1.0) ){
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStepBase, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }
                        for(size_t ii=0;ii<m_SaveFreq;ii++){
                            size_t iStep = m_NumberOfStepDonePreviously + iStepBase + ii;
                            IntegrateStepEta(iStep);
                            /*if( iStep % m_StatisticsFreq == 0 )*/ UpdateStatistics(iStep);
                        }
                        Record(m_NumberOfStepDonePreviously + iStepBase + m_SaveFreq - 1);

                        if( bTimer ){
                            if( igroup == ngroup - 1 || loopTimer.Tick(1.0) ){
                                loopTimer.Update(iStepBase + m_SaveFreq - 1, m_NumberOfStep + 1);
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStepBase + m_SaveFreq, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }
                    }
                }
            }else{
                if( m_Tau0 == 0.0 ){
                    LoopTimer loopTimer;
                    if( bTimer ){
                        loopTimer.Start();
                        fprintf(stderr, "\n");
                    }
                    for(size_t iStep0=1;iStep0<=m_NumberOfStep;iStep0++){
                        if( bTimer ){
                            if( iStep0 == 1 || loopTimer.Tick(1.0) ){
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStep0, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }

                        size_t iStep = m_NumberOfStepDonePreviously + iStep0;

                        IntegrateStep(iStep);
                        /*if( iStep % m_StatisticsFreq == 0 )*/ UpdateStatistics(iStep);
                        if( iStep % m_SaveFreq == 0 ) Record(iStep);

                        if( bTimer ){
                            if( iStep0 == m_NumberOfStep || loopTimer.Tick(1.0) ){
                                loopTimer.Update(iStep0, m_NumberOfStep + 1);
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStep0 + 1, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }
                    }

                }else{
                    LoopTimer loopTimer;
                    if( bTimer ){
                        loopTimer.Start();
                        fprintf(stderr, "\n");
                    }
                    for(size_t iStep0=1;iStep0<=m_NumberOfStep;iStep0++){
                        if( bTimer ){
                            if( iStep0 == 1 || loopTimer.Tick(1.0) ){
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStep0, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }

                        size_t iStep = m_NumberOfStepDonePreviously + iStep0;
                        IntegrateStepEta(iStep);
                        /*if( iStep % m_StatisticsFreq == 0 )*/ UpdateStatistics(iStep);
                        if( iStep % m_SaveFreq == 0 ) Record(iStep);

                        if( bTimer ){
                            if( iStep0 == m_NumberOfStep || loopTimer.Tick(1.0) ){
                                loopTimer.Update(iStep0, m_NumberOfStep + 1);
                                fprintf(stderr, "\33[F%10zd/%10zd %s\n", iStep0 + 1, m_NumberOfStep + 1, loopTimer.GetPChar());
                            }
                        }
                    }
                }
            }

            if( m_WriteRestartFile ){
                if( !WriteRestartFile(iTrial) ){
                    fprintf(stderr, "Error: can't open restart file on rank %10d\n", m_Rank);
                }
            }
            return true;
        }
        virtual void IntegrateInitializeVariables(){
            m_M = _mmSIMD_set1_pd(m_M0);
            m_T = _mmSIMD_set1_pd(m_T0);
            m_kB = _mmSIMD_set1_pd(m_kB0);
            m_Gamma = _mmSIMD_set1_pd(m_Gamma0);
            m_TimeStep = _mmSIMD_set1_pd(m_TimeStep0);

            m_MInverse = _mmSIMD_div_pd(m_One, m_M);
            m_MHalf = _mmSIMD_mul_pd(m_M, m_Half);
            m_HalfTimeStep = _mmSIMD_mul_pd(m_TimeStep, m_Half);
            m_kBInverse = _mmSIMD_div_pd(m_One, m_kB);

            double expMinusGammaTimeStep = std::exp(-m_Gamma0 * m_TimeStep0);
            m_ExpMinusGammaTimeStep = _mmSIMD_set1_pd(expMinusGammaTimeStep);

            double varianceOfRandomForceIntegration = 
                sqrt(m_kB0 * m_T0 / m_M0 * (1.0 - std::exp(-2.0 * m_Gamma0 * m_TimeStep0)));
            m_VarianceOfRandomForceIntegration = _mmSIMD_set1_pd(varianceOfRandomForceIntegration);  

            m_PeriodicLength = _mmSIMD_set1_pd(m_PeriodicLength0);
            // m_BarrierHeight = _mmSIMD_set1_pd(m_BarrierHeight0);

            if( m_bEqualProbability ){
                // equal probability
                m_MinusTwoTimeStepTauInverse0 = -2.0 * m_TimeStep0 / m_Tau0;
                double transitionProbability0 = 0.0;
                m_Tau = _mmSIMD_set1_pd(m_Tau0);
                if( m_Tau0 == 0.0 ){
                    transitionProbability0 = 0.0;
                }else{
                    // transitionProbability0 = 1.0 - exp(-m_TimeStep0 / m_Tau0);
                    transitionProbability0 = 0.5 * (1.0 - exp(-2.0 * m_TimeStep0 / m_Tau0)); // probability to the other state
                }
                m_TauInverse0 = 1.0 / m_Tau0;
                m_TauInverse = _mmSIMD_set1_pd(m_TauInverse0);

                m_TransitionProbability = _mmSIMD_set1_pd(transitionProbability0);

                m_ProbabilityOfFirstState0 = 0.5;

                m_ProbabilityOfFirstState = _mmSIMD_set1_pd(m_ProbabilityOfFirstState0);
            }else{
                // probability as energy
                m_MinusTwoTimeStepTauInverse0 = -(m_EquilibriumConstant + 1.0) * m_TimeStep0 / m_Tau0;
                double transitionProbabilityAsymm0[2] = {0.0, 0.0};
                m_Tau = _mmSIMD_set1_pd(m_Tau0); // k21
                if( m_Tau0 == 0.0 ){
                    transitionProbabilityAsymm0[0] = 0.0;
                    transitionProbabilityAsymm0[1] = 0.0;
                }else{
                    // transitionProbability0 = 1.0 - exp(-m_TimeStep0 / m_Tau0);
                    transitionProbabilityAsymm0[0] = 
                        m_EquilibriumConstant / (1.0 + m_EquilibriumConstant) * (1.0 - exp(m_MinusTwoTimeStepTauInverse0)); // probability from state 0 to state 1
                    transitionProbabilityAsymm0[1] = 
                        1.0 / (1.0 + m_EquilibriumConstant) * (1.0 - exp(m_MinusTwoTimeStepTauInverse0)); // probability from state 1 to state 0
                }
                m_TauInverse0 = 1.0 / m_Tau0;
                m_TauInverse = _mmSIMD_set1_pd(m_TauInverse0);

                m_TransitionProbabilityAsymm[0] = _mmSIMD_set1_pd(transitionProbabilityAsymm0[0]);
                m_TransitionProbabilityAsymm[1] = _mmSIMD_set1_pd(transitionProbabilityAsymm0[1]);

                m_ProbabilityOfFirstState0 = 1.0 / (1.0 + m_EquilibriumConstant);
                m_ProbabilityOfFirstState = _mmSIMD_set1_pd(m_ProbabilityOfFirstState0);
            }

            m_StateBarrierHeight_0 = _mmSIMD_set1_pd(m_BarrierHeight0[0]);
            m_StateBarrierHeight_1 = _mmSIMD_set1_pd(m_BarrierHeight0[1]);

            m_StateBarrierSlope_0 = _mmSIMD_set1_pd(m_BarrierHeight0[0] / m_SlopeHalfLength0);
            m_StateBarrierSlope_1 = _mmSIMD_set1_pd(m_BarrierHeight0[1] / m_SlopeHalfLength0);

            m_SlopeHalfLength = _mmSIMD_set1_pd(m_SlopeHalfLength0);

            m_SlopeBoundLow = _mmSIMD_set1_pd(-m_SlopeHalfLength0);
            m_SlopeBoundHigh = _mmSIMD_set1_pd(m_SlopeHalfLength0);

            m_HalfPeriodicLength = _mmSIMD_mul_pd(m_PeriodicLength, m_Half);
            m_PeriodicLengthInverse = _mmSIMD_div_pd(m_One, m_PeriodicLength);

            m_PeriodicLengthInverse0 = 1.0 / m_PeriodicLength0;

            m_HalfPeriodicLength0 = m_PeriodicLength0 * 0.5;

            // eta
            // if( m_Tau0 != 0.0 ){
            //     m_Q = _mmSIMD_set1_pd(m_Q0);
            //     m_Tau = _mmSIMD_set1_pd(m_Tau0);

            //     double ornsteinUhlenbeckVariance = sqrt(m_Q0 * (1.0 - std::exp(-2.0 / m_Tau0 * m_TimeStep0)));
            //     double ornsteinUnlenbackExpMinusGammaTimeStep = std::exp(-m_TimeStep0 / m_Tau0);

            //     m_OrensteinUhlenbeckVariance = _mmSIMD_set1_pd(ornsteinUhlenbeckVariance);
            //     m_OrensteinUnlenbackExpMinusGammaTimeStep = _mmSIMD_set1_pd(ornsteinUnlenbackExpMinusGammaTimeStep);

            // }


        };
        virtual void IntegrateStep0(){
            // x
            m_TrajectoryData->m_X = _mmSIMD_set1_pd(m_X0);
            InitializePeriodicIndices();

            // barrier state
            InitializeBarrierState();
            UpdateBarrierStateEtaParamter();

            // v
            if( m_RandomizeVelocity ){
                double velocitySigma = sqrt(m_kB0 * m_T0 / m_M0);
                __mSIMDd _velocitySigma = _mmSIMD_set1_pd(velocitySigma);
                m_TrajectoryData->m_V = _mmSIMD_mul_pd(_velocitySigma, m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe());
            }else{
                m_TrajectoryData->m_V = _mmSIMD_set1_pd(m_V0);
            }
            _mmSIMD_store_pd(m_TrajectoryData->m_V0Instance, m_TrajectoryData->m_V);

            // eta
            // if( m_Tau0 != 0.0 && m_RandomizeEta ){
            //     __mSIMDd sigma = _mmSIMD_set1_pd(sqrt(m_Q0));
            //     m_TrajectoryData->m_Eta = _mmSIMD_mul_pd(m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe(), sigma);
            // }else{
            //     m_TrajectoryData->m_Eta = _mmSIMD_set1_pd(0.0);
            // }
            // _mmSIMD_store_pd(m_TrajectoryData->m_Eta0Instance, m_TrajectoryData->m_Eta);

            // update m_E and _F
            EnergyAndForce();

            // initialize detail calculation
            InitializeStatistics();
        }

        virtual void IntegrateStep(size_t iStep){

            // v_n+1/3 = v_n + 1/m * F(qn) * dt/2
            // m_V = m_V + m_MInverse * m_F * m_HalfTimeStep;
            m_TrajectoryData->m_V = _mmSIMD_add_pd(m_TrajectoryData->m_V, _mmSIMD_mul_pd(_mmSIMD_mul_pd(m_MInverse, m_TrajectoryData->m_F), m_HalfTimeStep));
            // q_n+1/2 = q_n + v_n+1/3 * dt/2
            // m_X = m_X + m_V * m_HalfTimeStep;
            m_TrajectoryData->m_X = _mmSIMD_add_pd(m_TrajectoryData->m_X, _mmSIMD_mul_pd(m_TrajectoryData->m_V, m_HalfTimeStep));
            // v_n+2/3 = exp(-gamma*dt) * v_n+1/3 + sqrt(kB * T * (1-exp(-2*gamma*dt)) / m) * wn
            m_TrajectoryData->m_V = _mmSIMD_add_pd(
                _mmSIMD_mul_pd(m_ExpMinusGammaTimeStep, m_TrajectoryData->m_V), 
                _mmSIMD_mul_pd(m_VarianceOfRandomForceIntegration, m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe()));
            // q_n+1 = q_n+1/2 + v_n+2/3 * dt/2
            m_TrajectoryData->m_X = _mmSIMD_add_pd(m_TrajectoryData->m_X, _mmSIMD_mul_pd(m_TrajectoryData->m_V, m_HalfTimeStep));
            // v_n+1 = v_n+2/3 + 1/m * F(q_n+1) * dt/2

            UpdatePeriodicIndices(iStep);
            // even without transition, need to update current barrier when going to another barrier
            UpdateBarrierState(iStep);
            UpdateBarrierStateEtaParamter();

            // update m_E and _F
            EnergyAndForce();
            m_TrajectoryData->m_V = _mmSIMD_add_pd(m_TrajectoryData->m_V, _mmSIMD_mul_pd(_mmSIMD_mul_pd(m_MInverse, m_TrajectoryData->m_F), m_HalfTimeStep));

        }
        virtual void IntegrateStepEta(size_t iStep){

            // v_n+1/3 = v_n + 1/m * F(qn) * dt/2
            // m_V = m_V + m_MInverse * m_F * m_HalfTimeStep;
            m_TrajectoryData->m_V = _mmSIMD_add_pd(m_TrajectoryData->m_V, _mmSIMD_mul_pd(_mmSIMD_mul_pd(m_MInverse, m_TrajectoryData->m_F), m_HalfTimeStep));
            // q_n+1/2 = q_n + v_n+1/3 * dt/2
            // m_X = m_X + m_V * m_HalfTimeStep;
            m_TrajectoryData->m_X = _mmSIMD_add_pd(m_TrajectoryData->m_X, _mmSIMD_mul_pd(m_TrajectoryData->m_V, m_HalfTimeStep));
            // v_n+2/3 = exp(-gamma*dt) * v_n+1/3 + sqrt(kB * T * (1-exp(-2*gamma*dt)) / m) * wn
            m_TrajectoryData->m_V = _mmSIMD_add_pd(
                _mmSIMD_mul_pd(m_ExpMinusGammaTimeStep, m_TrajectoryData->m_V), 
                _mmSIMD_mul_pd(m_VarianceOfRandomForceIntegration, m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe()));
            // q_n+1 = q_n+1/2 + v_n+2/3 * dt/2
            m_TrajectoryData->m_X = _mmSIMD_add_pd(m_TrajectoryData->m_X, _mmSIMD_mul_pd(m_TrajectoryData->m_V, m_HalfTimeStep));
            // v_n+1 = v_n+2/3 + 1/m * F(q_n+1) * dt/2

            UpdatePeriodicIndices(iStep);

            // have to update the barrier here
            UpdateBarrierStateEta(iStep);
            UpdateBarrierStateEtaParamter();

            // update m_E and _F
            EnergyAndForce();
            m_TrajectoryData->m_V = _mmSIMD_add_pd(m_TrajectoryData->m_V, _mmSIMD_mul_pd(_mmSIMD_mul_pd(m_MInverse, m_TrajectoryData->m_F), m_HalfTimeStep));

        }
        virtual void Record(size_t iStep){
            if( m_SaveX ){
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                _mmSIMD_store_pd(buffer, m_TrajectoryData->m_X);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajX[i].push_back(buffer[i]);
                }
            }
            if( m_SaveV ){
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                _mmSIMD_store_pd(buffer, m_TrajectoryData->m_V);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajV[i].push_back(buffer[i]);
                }
            }
            if( m_SaveF ){
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                _mmSIMD_store_pd(buffer, m_TrajectoryData->m_F);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajF[i].push_back(buffer[i]);
                }
            }
            if( m_SaveE ){
                // PE
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                _mmSIMD_store_pd(buffer, m_TrajectoryData->m_E);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajPE[i].push_back(buffer[i]);
                }
                // KE
                __mSIMDd _ke = KineticEnergy(m_TrajectoryData->m_V);
                _mmSIMD_store_pd(buffer, _ke);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajKE[i].push_back(buffer[i]);
                }
            }
            if( m_SaveEta ){
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                _mmSIMD_store_pd(buffer, m_TrajectoryData->m_Eta);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajEta[i].push_back(buffer[i]);
                }
            }
            if( m_SaveState ){
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                _mmSIMD_store_pd(buffer, m_TrajectoryData->m_CurrentState);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajState[i].push_back(buffer[i]);
                }
            }
            if( m_SaveBarrier ){
                alignas(__mSIMDd) double buffer[SIMDWIDTH];
                __mSIMDd _barrier = m_TrajectoryData->m_Barrier;
                _mmSIMD_store_pd(buffer, _barrier);
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TrajBarrier[i].push_back(buffer[i]);
                }
            }
        }
        virtual void PrepareRecord(){
            m_TrajectoryData->Clear();
            // prepare record
            // X
            int numberOfFrame = CalculateNumberOfFrame(m_NumberOfStep, m_SaveFreq);
            for(int i=0;i<SIMDWIDTH;i++){
                m_TrajectoryData->m_TrajX[i].clear();
                m_TrajectoryData->m_TrajV[i].clear();
                m_TrajectoryData->m_TrajF[i].clear();
                m_TrajectoryData->m_TrajPE[i].clear();
                m_TrajectoryData->m_TrajKE[i].clear();
                m_TrajectoryData->m_TrajEta[i].clear();
                m_TrajectoryData->m_TrajState[i].clear();
                m_TrajectoryData->m_TrajBarrier[i].clear();

                if( m_SaveX ) m_TrajectoryData->m_TrajX[i].reserve(numberOfFrame);
                // V
                if( m_SaveV ) m_TrajectoryData->m_TrajV[i].reserve(numberOfFrame);
                // F
                if( m_SaveF ) m_TrajectoryData->m_TrajF[i].reserve(numberOfFrame);
                // E
                if( m_SaveE ){
                    m_TrajectoryData->m_TrajPE[i].reserve(numberOfFrame);
                    m_TrajectoryData->m_TrajKE[i].reserve(numberOfFrame);
                }
                // Eta
                if( m_SaveEta ){
                    m_TrajectoryData->m_TrajEta[i].reserve(numberOfFrame);
                }

                // barrier
                if( m_SaveBarrier ){
                    m_TrajectoryData->m_TrajBarrier[i].reserve(numberOfFrame);
                }
                // state
                if( m_SaveState ){
                    m_TrajectoryData->m_TrajState[i].reserve(numberOfFrame);
                }

                // for waiting time
                if( m_bWaitingTimeDistribution ){
                    m_TrajectoryData->m_BarrierTimeTick[i].reserve(100);
                }
            }
        }
        static size_t CalculateNumberOfFrame(size_t numberOfFrame, size_t freq){
            return numberOfFrame / freq + 1;
        }
        double CalculateNumberOfTrajectorySet(){
            double numberOfTrajectorySet = 0;
            if( m_SaveX ) numberOfTrajectorySet++;
            if( m_SaveV ) numberOfTrajectorySet++;
            if( m_SaveF ) numberOfTrajectorySet++;
            if( m_SaveE ) numberOfTrajectorySet += 2;
            if( m_SaveEta ) numberOfTrajectorySet++;
            if( m_SaveState ) numberOfTrajectorySet++;
            if( m_SaveBarrier ) numberOfTrajectorySet++;
            return numberOfTrajectorySet;
        }
        double CalculateMemoryUsage(){
            double size0 = CalculateNumberOfFrame(m_NumberOfStep, m_SaveFreq);
            return CalculateNumberOfTrajectorySet() * size0 * sizeof(double) * SIMDWIDTH;
        }
        int GetGaussianRandomNumberSIMDBankSize(){
            return m_GaussianRandomNumberSequenceSIMD.GetCoreBankSize();
        }
        int GetUniformRandomNumberSIMDBankSize(){
            return m_UniformRandomNumberSequenceSIMD.GetCoreBankSize();
        }
        int GetUniformRandomNumberBankSize(){
            return m_UniformRandomNumberSequence.GetCoreBankSize();
        }
        void InitializePeriodicIndices(){
            // update 
            __mSIMDd _reducedX = _mmSIMD_mul_pd(m_TrajectoryData->m_X, m_PeriodicLengthInverse);
            m_TrajectoryData->m_WellIndex = _mmSIMD_round_pd(_reducedX, 0x09);
            m_TrajectoryData->m_PeriodicIndex = _mmSIMD_round_pd(_reducedX, _MM_FROUND_TO_NEAREST_INT);

            m_TrajectoryData->m_WellIndexPrevious = m_TrajectoryData->m_WellIndex;
            m_TrajectoryData->m_PeriodicIndexPrevious = m_TrajectoryData->m_PeriodicIndex;
            m_TrajectoryData->m_bOnEnteringDifferentWell = 0;
            m_TrajectoryData->m_bOnEnteringDifferentPeriodicity = 0;

            // update stored indices
            _mmSIMD_store_pd(m_TrajectoryData->m_WellIndexStored, m_TrajectoryData->m_WellIndex);
            _mmSIMD_store_pd(m_TrajectoryData->m_PeriodicIndexStored, m_TrajectoryData->m_PeriodicIndex);
            _mmask8_store_SIMD(m_TrajectoryData->m_bOnEnteringDifferentWellStored, m_TrajectoryData->m_bOnEnteringDifferentWell);
            _mmask8_store_SIMD(m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored, m_TrajectoryData->m_bOnEnteringDifferentPeriodicity);

            memcpy(m_TrajectoryData->m_WellIndexPreviousStored, m_TrajectoryData->m_WellIndexStored, sizeof(double) * SIMDWIDTH);
            memcpy(m_TrajectoryData->m_PeriodicIndexPreviousStored, m_TrajectoryData->m_PeriodicIndexStored, sizeof(double) * SIMDWIDTH);

// printf("%10.0f aaa\n", m_TrajectoryData->m_PeriodicIndexStored[0]);
        }

        void UpdatePeriodicIndices(size_t iStep){
            // record previous 
            m_TrajectoryData->m_WellIndexPrevious = m_TrajectoryData->m_WellIndex;
            m_TrajectoryData->m_PeriodicIndexPrevious = m_TrajectoryData->m_PeriodicIndex;
            memcpy(m_TrajectoryData->m_WellIndexPreviousStored, m_TrajectoryData->m_WellIndexStored, sizeof(double) * SIMDWIDTH);
            memcpy(m_TrajectoryData->m_PeriodicIndexPreviousStored, m_TrajectoryData->m_PeriodicIndexStored, sizeof(double) * SIMDWIDTH);

            // update 
            __mSIMDd _reducedX = _mmSIMD_mul_pd(m_TrajectoryData->m_X, m_PeriodicLengthInverse);
            m_TrajectoryData->m_WellIndex = _mmSIMD_round_pd(_reducedX, 0x09);
            m_TrajectoryData->m_PeriodicIndex = _mmSIMD_round_pd(_reducedX, _MM_FROUND_TO_NEAREST_INT);

            // on entering event
            m_TrajectoryData->m_bOnEnteringDifferentWell = _mmSIMD_cmp_pd_mask(
                m_TrajectoryData->m_WellIndexPrevious, m_TrajectoryData->m_WellIndex, _CMP_NEQ_OQ);
            m_TrajectoryData->m_bOnEnteringDifferentPeriodicity = _mmSIMD_cmp_pd_mask( 
                m_TrajectoryData->m_PeriodicIndexPrevious, m_TrajectoryData->m_PeriodicIndex, _CMP_NEQ_OQ);

            // update stored indices
            _mmSIMD_store_pd(m_TrajectoryData->m_WellIndexStored, m_TrajectoryData->m_WellIndex);
            _mmSIMD_store_pd(m_TrajectoryData->m_PeriodicIndexStored, m_TrajectoryData->m_PeriodicIndex);
            _mmask8_store_SIMD(m_TrajectoryData->m_bOnEnteringDifferentWellStored, m_TrajectoryData->m_bOnEnteringDifferentWell);
            _mmask8_store_SIMD(m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored, m_TrajectoryData->m_bOnEnteringDifferentPeriodicity);
        }
        void InitializeBarrierState(){            

            // randomize initial state            
            if( m_RandomizeState0 ){
                __mSIMDd _randomNumber = m_UniformRandomNumberSequenceSIMD.GetRandomNumberSafe();       
                m_TrajectoryData->m_CurrentState = 
                    _mmSIMD_mask_blend_pd(
                        _mmSIMD_cmp_pd_mask(_randomNumber, m_ProbabilityOfFirstState, 2), m_One, m_Zero); // 2 : _CMP_LE_OS
            }else{
                m_TrajectoryData->m_CurrentState = _mmSIMD_set1_pd(m_InitialState0);
            }

            _mmSIMD_store_pd(m_TrajectoryData->m_BarrierState0Instance, m_TrajectoryData->m_CurrentState);
            m_TrajectoryData->m_PreviousState = m_TrajectoryData->m_CurrentState;

            const int leftmost = -100;
            const int rightmost = 100;
            m_TrajectoryData->m_BarrierStatus[leftmost];
            m_TrajectoryData->m_BarrierStatus[rightmost];

            // override initial state to the initial barrier
            for(int i=0;i<SIMDWIDTH;i++){
                int index = m_TrajectoryData->m_PeriodicIndexStored[i];

                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[index]);
                barrierStatus->m_State[i] = m_TrajectoryData->m_BarrierState0Instance[i];
            }

            m_TrajectoryData->m_StateCount[0] = _mmSIMD_set1_pd(0.0);
            m_TrajectoryData->m_StateCount[1] = _mmSIMD_set1_pd(0.0);

            // debug
            // for(int j=0;j<SIMDWIDTH;j++){
            //     char buffer[256];
            //     sprintf(buffer, "debug_initialBarrierStatus_%010d_%010d.log", m_Rank, j);
            //     FILE* fout = fopen(buffer, "w+");
            //     double average = 0;
            //     double count = 0;
            //     for(int i=leftmost;i<=rightmost;i++){
            //         BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[i]);
            //         fprintf(fout, "%10d %10d\n", i, barrierStatus->m_State[j]);
            //         average += barrierStatus->m_State[j];
            //         count++;
            //     }
            //     fprintf(fout, "# %13.7f %13.0f\n", average / count, count);
            //     fclose(fout);
            // }
        }
        void UpdateBarrierState(size_t iStep){
#ifndef SYNCHRONIZED
            // record previous 
            // m_TrajectoryData->m_PreviousState = m_TrajectoryData->m_CurrentState;

            // alignas(SIMDBYTESIZE) double previousState[SIMDWIDTH];
            // _mmSIMD_store_pd(previousState, m_TrajectoryData->m_PreviousState);

            // update barrier no matter to change barrier or not
            // due to complication of branching
            // not to update current state as it is static
            // __mSIMDd _uniformRandomNumber = m_UniformRandomNumberSequenceSIMD.GetRandomNumberSafe();
            // __mmask8 _toChange = 
            //     _mmSIMD_cmp_pd_mask(_uniformRandomNumber, m_TransitionProbability, 1); // prob < m_TransitionProbability
            // __mmask8 _condition = _mmSIMD_cmp_pd_mask(m_TrajectoryData->m_CurrentState, m_Zero, _CMP_EQ_OQ);
            // __mSIMDd _otherValue = _mmSIMD_mask_blend_pd(_condition, m_Zero, m_One);
            // m_TrajectoryData->m_CurrentState = _mmSIMD_mask_blend_pd(_toChange, m_TrajectoryData->m_CurrentState, _otherValue);

            // alignas(SIMDBYTESIZE) double uniformRandomNumber[SIMDWIDTH];
            alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];

            // _mmSIMD_store_pd(uniformRandomNumber, _uniformRandomNumber);
            _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);

            // determine barrier index
            for(int i=0;i<SIMDWIDTH;i++){
                if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                    // save status of the previous barrier
                    int indexPrevious = m_TrajectoryData->m_PeriodicIndexPreviousStored[i];
                    BarrierStatus* barrierStatusPrevious = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[indexPrevious]);

                    barrierStatusPrevious->m_State[i] = currentState[i];
                    // now advance this barrier at this timestep so iStep
                    barrierStatusPrevious->m_TimeStep[i] = iStep; 
                    // imagine the same barrier is saved and taken at the same step,
                    // the saved time step is the previous one so iStep - 1
                    // barrierStatusPrevious->m_TimeStep[i] = iStep - 1; 

                    // take status of the current barrier
                    int indexCurrent = m_TrajectoryData->m_PeriodicIndexStored[i];
                    BarrierStatus* barrierStatusCurrent = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[indexCurrent]);
                    // previousState[i] = barrierStatusCurrent->m_State[i];
                    currentState[i] = barrierStatusCurrent->m_State[i];

                    // if( i == 7 ) printf("%13zd %13zd %2d(%d)->%2d(%d) \n", 
                    //     barrierStatusCurrent->m_TimeStep[i], iStep, 
                    //     indexPrevious, barrierStatusPrevious->m_State[i], indexCurrent, (int)currentState[i]);
                }

                // debug
                // if( i==0 && m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                //     printf("%10zd %3d -> %3d %13.7f %13.7f %13.7f %13.7f %1d %1.0f diff\n", 
                //         iStep, indexPrevious, indexCurrent, m_TrajectoryData->m_X[0], timeLap, transitionProbability, m_TransitionProbability[0], 
                //         barrierStatusPrevious->m_State[i], currentState[i]);
                // }else if( i == 0 ){
                //     printf("%10zd %3d -> %3d %13.7f %13.7f %13.7f %13.7f %1d %1.0f\n", 
                //         iStep, indexPrevious, indexCurrent, m_TrajectoryData->m_X[0], timeLap, transitionProbability, m_TransitionProbability[0],
                //         barrierStatusPrevious->m_State[i], currentState[i]);
                // }
            }
            // update m_CurrentState using changed barrier index
            m_TrajectoryData->m_CurrentState = _mmSIMD_load_pd(currentState);
            
#endif
        }
        void UpdateBarrierStateEta(size_t iStep){
            // record previous 
            // m_TrajectoryData->m_PreviousState = m_TrajectoryData->m_CurrentState;

            // alignas(SIMDBYTESIZE) double previousState[SIMDWIDTH];
            // _mmSIMD_store_pd(previousState, m_TrajectoryData->m_PreviousState);

            if( m_bEqualProbability ){
                // equal probability

                // update barrier no matter to change barrier or not
                // due to complication of branching
                __mSIMDd _uniformRandomNumber = m_UniformRandomNumberSequenceSIMD.GetRandomNumberSafe();
                __mmask8 _toChange = 
                    _mmSIMD_cmp_pd_mask(_uniformRandomNumber, m_TransitionProbability, 1); // prob < m_TransitionProbability
                __mmask8 _condition = _mmSIMD_cmp_pd_mask(m_TrajectoryData->m_CurrentState, m_Zero, _CMP_EQ_OQ);
                __mSIMDd _otherValue = _mmSIMD_mask_blend_pd(_condition, m_Zero, m_One);
                m_TrajectoryData->m_CurrentState = _mmSIMD_mask_blend_pd(_toChange, m_TrajectoryData->m_CurrentState, _otherValue);

#ifndef SYNCHRONIZED
                // alignas(SIMDBYTESIZE) double uniformRandomNumber[SIMDWIDTH];
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];

                // _mmSIMD_store_pd(uniformRandomNumber, _uniformRandomNumber);
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);

                // determine barrier index
                for(int i=0;i<SIMDWIDTH;i++){
                    if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        // save status of the previous barrier
                        int indexPrevious = m_TrajectoryData->m_PeriodicIndexPreviousStored[i];
                        BarrierStatus* barrierStatusPrevious = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[indexPrevious]);

                        barrierStatusPrevious->m_State[i] = currentState[i];
                        // now advance this barrier at this timestep so iStep
                        barrierStatusPrevious->m_TimeStep[i] = iStep; 
                        // imagine the same barrier is saved and taken at the same step,
                        // the saved time step is the previous one so iStep - 1
                        // barrierStatusPrevious->m_TimeStep[i] = iStep - 1; 

                        // take status of the current barrier
                        int indexCurrent = m_TrajectoryData->m_PeriodicIndexStored[i];
                        BarrierStatus* barrierStatusCurrent = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[indexCurrent]);
                        // previousState[i] = barrierStatusCurrent->m_State[i];
                        double previousState = barrierStatusCurrent->m_State[i];
                        double timeLap = iStep - barrierStatusCurrent->m_TimeStep[i];
                        double transitionProbability = 0.5 * (1.0 - exp(m_MinusTwoTimeStepTauInverse0 * timeLap)); // probability to the other state
                        double uniformRandomNumber = m_UniformRandomNumberSequence.GetRandomNumber();
                        if( uniformRandomNumber < transitionProbability ){
                            if( previousState == 1.0 ) currentState[i] = 0.0;
                            else currentState[i] = 1.0;
                        }else{
                            currentState[i] = previousState;
                        }
                        // if( i == 7 ) printf("%13zd %13zd %13.7f %13.7f %2d(%d)->%2d(%d) %13.7f %d (%13.7f %13.7f)\n", 
                        //     barrierStatusCurrent->m_TimeStep[i], iStep, timeLap, transitionProbability, 
                        //     indexPrevious, barrierStatusPrevious->m_State[i], indexCurrent, (int)currentState[i], previousState, previousState == 1.0,
                        //     uniformRandomNumber, transitionProbability);
                    }

                    // debug
                    // if( i==0 && m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                    //     printf("%10zd %3d -> %3d %13.7f %13.7f %13.7f %13.7f %1d %1.0f diff\n", 
                    //         iStep, indexPrevious, indexCurrent, m_TrajectoryData->m_X[0], timeLap, transitionProbability, m_TransitionProbability[0], 
                    //         barrierStatusPrevious->m_State[i], currentState[i]);
                    // }else if( i == 0 ){
                    //     printf("%10zd %3d -> %3d %13.7f %13.7f %13.7f %13.7f %1d %1.0f\n", 
                    //         iStep, indexPrevious, indexCurrent, m_TrajectoryData->m_X[0], timeLap, transitionProbability, m_TransitionProbability[0],
                    //         barrierStatusPrevious->m_State[i], currentState[i]);
                    // }
                }

                // update m_CurrentState using changed barrier index
                m_TrajectoryData->m_CurrentState = _mmSIMD_load_pd(currentState);

#endif // SYNCHRONIZED
            }else{
                // probability as energy
                // update barrier no matter to change barrier or not
                // due to complication of branching
                __mmask8 _condition = _mmSIMD_cmp_pd_mask(m_TrajectoryData->m_CurrentState, m_One, _CMP_EQ_OQ);
                __mSIMDd transitionProbability = _mmSIMD_mask_blend_pd(_condition, m_TransitionProbabilityAsymm[0], m_TransitionProbabilityAsymm[1]);
                __mSIMDd _uniformRandomNumber = m_UniformRandomNumberSequenceSIMD.GetRandomNumberSafe();
                __mmask8 _toChange = 
                    _mmSIMD_cmp_pd_mask(_uniformRandomNumber, transitionProbability, 1); // prob < m_TransitionProbability
                // __mmask8 _condition = _mmSIMD_cmp_pd_mask(m_TrajectoryData->m_CurrentState, m_Zero, _CMP_EQ_OQ);
                // __mSIMDd _otherValue = _mmSIMD_mask_blend_pd(_condition, m_Zero, m_One);
                __mSIMDd _otherValue = _mmSIMD_mask_blend_pd(_condition, m_One, m_Zero);
                m_TrajectoryData->m_CurrentState = _mmSIMD_mask_blend_pd(_toChange, m_TrajectoryData->m_CurrentState, _otherValue);

#ifndef SYNCHRONIZED
                // alignas(SIMDBYTESIZE) double uniformRandomNumber[SIMDWIDTH];
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];

                // _mmSIMD_store_pd(uniformRandomNumber, _uniformRandomNumber);
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);

                // determine barrier index
                for(int i=0;i<SIMDWIDTH;i++){
                    if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        // save status of the previous barrier
                        int indexPrevious = m_TrajectoryData->m_PeriodicIndexPreviousStored[i];
                        BarrierStatus* barrierStatusPrevious = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[indexPrevious]);

                        barrierStatusPrevious->m_State[i] = currentState[i];
                        // now advance this barrier at this timestep so iStep
                        barrierStatusPrevious->m_TimeStep[i] = iStep; 
                        // imagine the same barrier is saved and taken at the same step,
                        // the saved time step is the previous one so iStep - 1
                        // barrierStatusPrevious->m_TimeStep[i] = iStep - 1; 

                        // take status of the current barrier
                        int indexCurrent = m_TrajectoryData->m_PeriodicIndexStored[i];
                        BarrierStatus* barrierStatusCurrent = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[indexCurrent]);
                        // previousState[i] = barrierStatusCurrent->m_State[i];
                        double previousState = barrierStatusCurrent->m_State[i];
                        double timeLap = iStep - barrierStatusCurrent->m_TimeStep[i];
                        double transitionProbability;
                        if( previousState == 0.0 ){
                            transitionProbability = m_EquilibriumConstant / (1.0 + m_EquilibriumConstant) * (1.0 - exp(m_MinusTwoTimeStepTauInverse0 * timeLap));
                        }else{
                            transitionProbability = 1.0 / (1.0 + m_EquilibriumConstant) * (1.0 - exp(m_MinusTwoTimeStepTauInverse0 * timeLap));
                        }
                        // double transitionProbability = 0.5 * (1.0 - exp(m_MinusTwoTimeStepTauInverse0 * timeLap)); // probability to the other state
                        double uniformRandomNumber = m_UniformRandomNumberSequence.GetRandomNumber();
                        if( uniformRandomNumber < transitionProbability ){
                            if( previousState == 1.0 ) currentState[i] = 0.0;
                            else currentState[i] = 1.0;
                        }else{
                            currentState[i] = previousState;
                        }
                        // if( i == 7 ) printf("%13zd %13zd %13.7f %13.7f %2d(%d)->%2d(%d) %13.7f %d (%13.7f %13.7f)\n", 
                        //     barrierStatusCurrent->m_TimeStep[i], iStep, timeLap, transitionProbability, 
                        //     indexPrevious, barrierStatusPrevious->m_State[i], indexCurrent, (int)currentState[i], previousState, previousState == 1.0,
                        //     uniformRandomNumber, transitionProbability);
                    }

                    // debug
                    // if( i==0 && m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                    //     printf("%10zd %3d -> %3d %13.7f %13.7f %13.7f %13.7f %1d %1.0f diff\n", 
                    //         iStep, indexPrevious, indexCurrent, m_TrajectoryData->m_X[0], timeLap, transitionProbability, m_TransitionProbability[0], 
                    //         barrierStatusPrevious->m_State[i], currentState[i]);
                    // }else if( i == 0 ){
                    //     printf("%10zd %3d -> %3d %13.7f %13.7f %13.7f %13.7f %1d %1.0f\n", 
                    //         iStep, indexPrevious, indexCurrent, m_TrajectoryData->m_X[0], timeLap, transitionProbability, m_TransitionProbability[0],
                    //         barrierStatusPrevious->m_State[i], currentState[i]);
                    // }
                }

                // update m_CurrentState using changed barrier index
                m_TrajectoryData->m_CurrentState = _mmSIMD_load_pd(currentState);

#endif // SYNCHRONIZED

            } // end if equal probability
        }
        virtual void UpdateBarrierStateEtaParamter(){
            __mmask8 _condition = _mmSIMD_cmp_pd_mask(m_TrajectoryData->m_CurrentState, m_One, _CMP_EQ_OQ);
            m_TrajectoryData->m_CurrentBarrierHeight = _mmSIMD_mask_blend_pd(_condition, m_StateBarrierHeight_0, m_StateBarrierHeight_1);
            // multiply by -1
            m_TrajectoryData->m_SlopePositive = _mmSIMD_mask_blend_pd(_condition, m_StateBarrierSlope_0, m_StateBarrierSlope_1);
            m_TrajectoryData->m_SlopeNegative = _mmSIMD_sub_pd(m_Zero, m_TrajectoryData->m_SlopePositive);

            // cummulate state count
            //                                                                    false=0,   true=1
            m_TrajectoryData->m_StateCount[0] = _mmSIMD_mask_blend_pd(_condition, _mmSIMD_add_pd(m_TrajectoryData->m_StateCount[0], m_One), m_TrajectoryData->m_StateCount[0]);
            m_TrajectoryData->m_StateCount[1] = _mmSIMD_mask_blend_pd(_condition, m_TrajectoryData->m_StateCount[1], _mmSIMD_add_pd(m_TrajectoryData->m_StateCount[1], m_One));
        }
        void InitializeStatistics(){
            // this routine is not called when restarted

            if( m_bWaitingTimeDistribution ){
                for(int i=0;i<SIMDWIDTH;i++){
                    m_TrajectoryData->m_TimeTickPrevious[i] = 0;
                    m_TrajectoryData->m_BarrierPassed[i] = false;
                }
            }

            if( m_bFirstPassageTime ){
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);
                
                for(int i=0;i<SIMDWIDTH;i++){
                    // from bottom to top
                    m_TrajectoryData->m_FirstPassageTime[i].clear();

                    if( Math::Modulo(m_X0, m_PeriodicLength0) == m_HalfPeriodicLength0 ){
                        m_TrajectoryData->m_FirstPassageTimeLastTick[i] = 0.0;
                        m_TrajectoryData->m_FirstPassageTimeToUpdate[i] = false;
                        if( m_TrajectoryData->m_WellIndexStored[i] == m_TrajectoryData->m_PeriodicIndexStored[i] ){
                            // on \ of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierState[i][0] = currentState[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] + 1]);
                                m_TrajectoryData->m_InitialBarrierState[i][1] = barrierStatus->m_State[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i] + 1;
                            }
                        }else{
                            // on / of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] - 1]);
                                m_TrajectoryData->m_InitialBarrierState[i][0] = barrierStatus->m_State[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i] - 1;
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierState[i][1] = currentState[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                        }
                        // debug
                        // if( i == 0 ){
                        //     fprintf(stderr, "%10d : %.0f , %10d : %.0f\n\n", 
                        //         m_TrajectoryData->m_InitialBarrierIndex[i][0], m_TrajectoryData->m_InitialBarrierState[i][0], 
                        //         m_TrajectoryData->m_InitialBarrierIndex[i][1], m_TrajectoryData->m_InitialBarrierState[i][1]);
                        // }
                    }else{
                        m_TrajectoryData->m_FirstPassageTimeLastTick[i] = 0.0;
                        m_TrajectoryData->m_FirstPassageTimeToUpdate[i] = true;
                        m_TrajectoryData->m_InitialBarrierState[i][0] = 0.0;
                    }

                    // from bottom to another bottom
                    m_TrajectoryData->m_FirstPassageTime2[i].clear();

                    if( Math::Modulo(m_X0, m_PeriodicLength0) == m_HalfPeriodicLength0 ){
                        m_TrajectoryData->m_FirstPassageTimeLastTick2[i] = 0.0;
                        m_TrajectoryData->m_FirstPassageTimeToUpdate2[i] = false;
                        if( m_TrajectoryData->m_WellIndexStored[i] == m_TrajectoryData->m_PeriodicIndexStored[i] ){
                            // on \ of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] + 1]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i] + 1;
                            }
                        }else{
                            // on / of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] - 1]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i] - 1;
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                        }
                        // debug
                        // if( i == 0 ){
                        //     fprintf(stderr, "%10d : %.0f , %10d : %.0f\n\n", 
                        //         m_TrajectoryData->m_InitialBarrierIndex2[i][0], m_TrajectoryData->m_InitialBarrierState[i][0], 
                        //         m_TrajectoryData->m_InitialBarrierIndex2[i][1], m_TrajectoryData->m_InitialBarrierState[i][1]);
                        // }
                    }else{
                        m_TrajectoryData->m_FirstPassageTimeLastTick2[i] = 0.0;
                        m_TrajectoryData->m_FirstPassageTimeToUpdate2[i] = true;
                    }
                }

            }

            if( m_bTransitionPathTime ){
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);
                
                for(int i=0;i<SIMDWIDTH;i++){
                    // from bottom to top
                    m_TrajectoryData->m_TransitionPathTime[i].clear();

                    if( Math::Modulo(m_X0, m_PeriodicLength0) == m_HalfPeriodicLength0 ){
                        m_TrajectoryData->m_TransitionPathTimeLastTick[i] = 0.0;
                        m_TrajectoryData->m_TransitionPathTimeToUpdate[i] = false;
                    }else{
                        m_TrajectoryData->m_TransitionPathTimeLastTick[i] = 0.0;
                        m_TrajectoryData->m_TransitionPathTimeToUpdate[i] = true;
                    }

                    // from bottom to another bottom
                    m_TrajectoryData->m_TransitionPathTime2[i].clear();

                    if( Math::Modulo(m_X0, m_PeriodicLength0) == m_HalfPeriodicLength0 ){
                        m_TrajectoryData->m_TransitionPathTimeLastTick2[i] = 0.0;
                        m_TrajectoryData->m_TransitionPathTimeToUpdate2[i] = false;
                        m_TrajectoryData->m_TransitionPathTimeInitialWell2[i] = m_TrajectoryData->m_WellIndexStored[i];
                    }else{
                        m_TrajectoryData->m_TransitionPathTimeLastTick2[i] = 0.0;
                        m_TrajectoryData->m_TransitionPathTimeToUpdate2[i] = true;
                        m_TrajectoryData->m_TransitionPathTimeInitialWell2[i] = m_TrajectoryData->m_WellIndexStored[i];
                    }
                }

            }

            if( m_bFirstPassageTimeCG ){
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);
                
                for(int i=0;i<SIMDWIDTH;i++){
                    for(int j=0;j<m_TrajectoryData->m_FirstPassageTimeCGSize.size();j++){
                        // // from bottom to top
                        // m_TrajectoryData->m_FirstPassageTimeCG[i].clear();

                        // if( Math::Modulo(m_X0, m_PeriodicLength0) == m_HalfPeriodicLength0 &&
                        //     (int)round(m_TrajectoryData->m_WellIndexStored[i]) % m_FirstPassageTimeCGSize == 0 ){
                        //     m_TrajectoryData->m_FirstPassageTimeLastTickCG[i] = 0.0;
                        //     m_TrajectoryData->m_FirstPassageTimeToUpdateCG[i] = false;
                        // }else{
                        //     m_TrajectoryData->m_FirstPassageTimeLastTickCG[i] = 0.0;
                        //     m_TrajectoryData->m_FirstPassageTimeToUpdateCG[i] = true;
                        // }

                        // from bottom to another bottom
                        m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i].clear();

                        if( Math::Modulo(m_X0, m_PeriodicLength0) == m_HalfPeriodicLength0 &&
                            (int)round(m_TrajectoryData->m_WellIndexStored[i]) % m_TrajectoryData->m_FirstPassageTimeCGSize[j] == 0 ){
                            m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i] = 0.0;
                            m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeToUpdateCG2[i] = false;
                            m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i] = m_TrajectoryData->m_WellIndexStored[i];
                        }else{
                            m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i] = 0.0;
                            m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeToUpdateCG2[i] = true;
                            m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i] = m_TrajectoryData->m_WellIndexStored[i];
                        }
                    }
                }

            }

        }
        void UpdateStatistics(size_t iStep){
            // update well and periodicity indices

            if( m_bWaitingTimeDistribution ){
                for(int i=0;i<SIMDWIDTH;i++){
                    if( m_TrajectoryData->m_bOnEnteringDifferentWellStored[i] ){
                        m_TrajectoryData->m_TimeTickPrevious[i] = iStep;
                        m_TrajectoryData->m_BarrierPassed[i] = true;
                    }else if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        if( m_TrajectoryData->m_BarrierPassed[i] ){
                            // zero means the particle has never crossed a well boundary, i.e. a barrier
                            m_TrajectoryData->m_BarrierTimeTick[i].push_back(m_TrajectoryData->m_TimeTickPrevious[i]);
                            m_TrajectoryData->m_BarrierPassed[i] = false;
                        }
                    }
                }
            }

            if( m_bFirstPassageTime ){
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);
                
                for(int i=0;i<SIMDWIDTH;i++){
                    // from bottom to top
                    if( m_TrajectoryData->m_FirstPassageTimeToUpdate[i] && 
                        m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        // hit starting line after barrier, so renew
                        m_TrajectoryData->m_FirstPassageTimeToUpdate[i] = false;
                        m_TrajectoryData->m_FirstPassageTimeLastTick[i] = m_TimeStep0 * iStep;

                        if( m_TrajectoryData->m_WellIndexStored[i] == m_TrajectoryData->m_PeriodicIndexStored[i] ){
                            // on \ of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierState[i][0] = currentState[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] + 1]);
                                m_TrajectoryData->m_InitialBarrierState[i][1] = barrierStatus->m_State[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i] + 1;
                            }
                        }else{
                            // on / of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] - 1]);
                                m_TrajectoryData->m_InitialBarrierState[i][0] = barrierStatus->m_State[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i] - 1;
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierState[i][1] = currentState[i];
                                m_TrajectoryData->m_InitialBarrierIndex[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                        }

                        // if( i==0){
                        //     fprintf(stdout, "update %23.16E %13.7f\n", m_TimeStep0 * iStep, ((double*)&m_TrajectoryData->m_X)[0]);
                        // }
                    }
                    if( !m_TrajectoryData->m_FirstPassageTimeToUpdate[i] && 
                        m_TrajectoryData->m_bOnEnteringDifferentWellStored[i] ){
                        // pass barrier
                        m_TrajectoryData->m_FirstPassageTimeToUpdate[i] = true;
                        m_TrajectoryData->m_FirstPassageTime[i].push_back(m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeLastTick[i]);

                        if( m_TrajectoryData->m_PeriodicIndexStored[i] == m_TrajectoryData->m_InitialBarrierIndex[i][0] ){
                            m_TrajectoryData->m_InitialBarrierList[i].push_back(m_TrajectoryData->m_InitialBarrierState[i][0]);
                        }else{
                            m_TrajectoryData->m_InitialBarrierList[i].push_back(m_TrajectoryData->m_InitialBarrierState[i][1]);
                        }
                        m_TrajectoryData->m_InstantBarrierList[i].push_back(currentState[i]);
                        // if( i==0){
                        //     fprintf(stdout, "pass %23.16E %13.7f %13.7f\n", 
                        //         m_TimeStep0 * iStep, ((double*)&m_TrajectoryData->m_X)[0], 
                        //         m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeLastTick[i]);
                        // }
                    }
                    // from bottom to another bottom
                    if( m_TrajectoryData->m_FirstPassageTimeToUpdate2[i] && 
                        m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        // hit starting line after barrier, so renew
                        m_TrajectoryData->m_FirstPassageTimeToUpdate2[i] = false;
                        m_TrajectoryData->m_FirstPassageTimeLastTick2[i] = m_TimeStep0 * iStep;

                        if( m_TrajectoryData->m_WellIndexStored[i] == m_TrajectoryData->m_PeriodicIndexStored[i] ){
                            // on \ of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] + 1]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i] + 1;
                            }
                        }else{
                            // on / of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] - 1]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i] - 1;
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                        }

                        // if( i==0){
                        //     fprintf(stdout, "update %23.16E %13.7f\n", m_TimeStep0 * iStep, ((double*)&m_TrajectoryData->m_X)[0]);
                        // }
                    }
                    if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] && 
                        !m_TrajectoryData->m_FirstPassageTimeToUpdate2[i] && 
                        ( m_TrajectoryData->m_InitialBarrierIndex2[i][0] != m_TrajectoryData->m_PeriodicIndexStored[i] &&
                          m_TrajectoryData->m_InitialBarrierIndex2[i][1] != m_TrajectoryData->m_PeriodicIndexStored[i] ) ){
                        // pass another bottom

                        // start counting immediately
                        // m_TrajectoryData->m_FirstPassageTimeToUpdate2[i] = true;
                        m_TrajectoryData->m_FirstPassageTime2[i].push_back(m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeLastTick2[i]);

                        // debug
                        // if( i==0){
                        //     fprintf(stdout, "pass %23.16E %13.7f %13.7f i %3d %3d n %3d\n", 
                        //         m_TimeStep0 * iStep, ((double*)&m_TrajectoryData->m_X)[0], 
                        //         m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeLastTick2[i],
                        //         m_TrajectoryData->m_InitialBarrierIndex2[i][0],
                        //         m_TrajectoryData->m_InitialBarrierIndex2[i][1],
                        //         (int)m_TrajectoryData->m_PeriodicIndexStored[i]);
                        // }

                        // restart
                        m_TrajectoryData->m_FirstPassageTimeLastTick2[i] = m_TimeStep0 * iStep;
                        if( m_TrajectoryData->m_WellIndexStored[i] == m_TrajectoryData->m_PeriodicIndexStored[i] ){
                            // on \ of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] + 1]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i] + 1;
                            }
                        }else{
                            // on / of \/
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i] - 1]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][0] = m_TrajectoryData->m_PeriodicIndexStored[i] - 1;
                            }
                            {
                                BarrierStatus* barrierStatus = static_cast<BarrierStatus*>(m_TrajectoryData->m_BarrierStatus[m_TrajectoryData->m_PeriodicIndexStored[i]]);
                                m_TrajectoryData->m_InitialBarrierIndex2[i][1] = m_TrajectoryData->m_PeriodicIndexStored[i];
                            }
                        }

                    }
                }
            }

            if( m_bTransitionPathTime ){
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);
                
                for(int i=0;i<SIMDWIDTH;i++){
                    // from bottom to top
                    if( /*m_TrajectoryData->m_TransitionPathTimeToUpdate[i] && */ // update anyway when hitting the bottom
                        m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        // hit starting line after barrier, so renew
                        m_TrajectoryData->m_TransitionPathTimeToUpdate[i] = false;
                        m_TrajectoryData->m_TransitionPathTimeLastTick[i] = m_TimeStep0 * iStep;
                    }
                    if( !m_TrajectoryData->m_TransitionPathTimeToUpdate[i] && 
                        m_TrajectoryData->m_bOnEnteringDifferentWellStored[i] ){
                        // pass barrier
                        m_TrajectoryData->m_TransitionPathTimeToUpdate[i] = true;
                        m_TrajectoryData->m_TransitionPathTime[i].push_back(m_TimeStep0 * iStep - m_TrajectoryData->m_TransitionPathTimeLastTick[i]);
                        // fprintf(stderr, "traversing1 %1d %23.16f %23.16f %23.16f\n\n", 
                        //     i, m_TimeStep0 * iStep, m_TrajectoryData->m_TransitionPathTimeLastTick[i],
                        //     m_TimeStep0 * iStep - m_TrajectoryData->m_TransitionPathTimeLastTick[i]);
                    }
                    
                    // from bottom to another bottom
                    if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        // hitting the bottom
                        if( m_TrajectoryData->m_TransitionPathTimeToUpdate2[i] ){
                            // not initialized
                            m_TrajectoryData->m_TransitionPathTimeToUpdate2[i] = false;
                        }else{
                            // initialized
                            if( m_TrajectoryData->m_WellIndexStored[i] != m_TrajectoryData->m_TransitionPathTimeInitialWell2[i] ){
                                // onto a different well
                                m_TrajectoryData->m_TransitionPathTime2[i].push_back(m_TimeStep0 * iStep - m_TrajectoryData->m_TransitionPathTimeLastTick2[i]);
                                // fprintf(stderr, "traversing2 %1d %23.16f %23.16f %23.16f\n\n", 
                                //     i, m_TimeStep0 * iStep, m_TrajectoryData->m_TransitionPathTimeLastTick2[i],
                                //     m_TimeStep0 * iStep - m_TrajectoryData->m_TransitionPathTimeLastTick2[i]);
                            }else{
                                // on the same well
                            }
                        }
                        m_TrajectoryData->m_TransitionPathTimeLastTick2[i] = m_TimeStep0 * iStep;
                        m_TrajectoryData->m_TransitionPathTimeInitialWell2[i] = m_TrajectoryData->m_WellIndexStored[i];
                    }
                }
            }

            if( m_bFirstPassageTimeCG ){
                alignas(SIMDBYTESIZE) double currentState[SIMDWIDTH];
                _mmSIMD_store_pd(currentState, m_TrajectoryData->m_CurrentState);
                
                for(int i=0;i<SIMDWIDTH;i++){
                    for(int j=0;j<m_TrajectoryData->m_FirstPassageTimeCGSize.size();j++){
                        // from bottom to top
                        // if( /*m_TrajectoryData->m_FirstPassageTimeCGToUpdate[i] && */ // update anyway when hitting the bottom
                        //     m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] ){
                        //     // hit starting line after barrier, so renew
                        //     m_TrajectoryData->m_FirstPassageTimeToUpdateCG[i] = false;
                        //     m_TrajectoryData->m_FirstPassageTimeLastTickCG[i] = m_TimeStep0 * iStep;
                        // }
                        // if( !m_TrajectoryData->m_FirstPassageTimeToUpdateCG[i] && 
                        //     m_TrajectoryData->m_bOnEnteringDifferentWellStored[i] ){
                        //     // pass barrier
                        //     m_TrajectoryData->m_FirstPassageTimeToUpdateCG[i] = true;
                        //     m_TrajectoryData->m_FirstPassageTimeCG[i].push_back(m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeLastTickCG[i]);
                        //     // fprintf(stderr, "traversing1 %1d %23.16f %23.16f %23.16f\n\n", 
                        //     //     i, m_TimeStep0 * iStep, m_TrajectoryData->m_FirstPassageTimeCGLastTick[i],
                        //     //     m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeCGLastTick[i]);
                        // }
                        
                        // from bottom to another bottom
                        if( m_TrajectoryData->m_bOnEnteringDifferentPeriodicityStored[i] &&
                            (int)round(m_TrajectoryData->m_WellIndexStored[i]) % m_TrajectoryData->m_FirstPassageTimeCGSize[j] == 0 ){
                            // hitting the bottom
                            if( m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeToUpdateCG2[i] ){
                                // not initialized
                                m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeToUpdateCG2[i] = false;
                                m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i] = m_TimeStep0 * iStep;
                                m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i] = m_TrajectoryData->m_WellIndexStored[i];
                            }else{
                                // initialized
                                if( m_TrajectoryData->m_WellIndexStored[i] != m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i] ){
                                    // onto a different well
                                    m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeCG2[i].push_back(
                                        m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i]);
                                    // fprintf(stderr, "fptcg2 %1d %23.16f %23.16f %23.16f %10d -> %10.0f\n\n", 
                                    //     i, m_TimeStep0 * iStep, m_TrajectoryData->m_FirstPassageTimeLastTickCG2[i],
                                    //     m_TimeStep0 * iStep - m_TrajectoryData->m_FirstPassageTimeLastTickCG2[i],
                                    //     m_TrajectoryData->m_FirstPassageTimeInitialWellCG2[i], m_TrajectoryData->m_WellIndexStored[i]);
                                    m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeLastTickCG2[i] = m_TimeStep0 * iStep;
                                    m_TrajectoryData->m_FirstPassageTimeCGData[j].m_FirstPassageTimeInitialWellCG2[i] = m_TrajectoryData->m_WellIndexStored[i];
                                }else{
                                    // on the same well
                                }
                            }
                        }
                    }
                }
            }

        }
    public:
        // virtual double GetWellNumber(double x){
        //     return floor((x/* - 0.0*/) * m_PeriodicLengthInverse0);
        // }
        // virtual double GetBarrierNumber(double x){
        //     return round((x/* - 0.0*/) * m_PeriodicLengthInverse0);
        // }
        // virtual double GetXInWell(double x){
        //     double xp = x - (floor((x/* - 0.0*/) * m_PeriodicLengthInverse0) + 0.5) * m_PeriodicLength0;
        //     return xp;
        // }
        // virtual double GetXP(double x){
        //     double xp = x - round((x/* - 0.0*/) * m_PeriodicLengthInverse0) * m_PeriodicLength0;
        //     return xp;
        // }
        // virtual double GetXPShift(double x){
        //     return - round((x/* - 0.0*/) * m_PeriodicLengthInverse0) * m_PeriodicLength0;
        // }
        // virtual bool OnBarrier(double x){
        //     double xp = x - round((x/* - 0.0*/) * m_PeriodicLengthInverse0) * m_PeriodicLength0;
        //     if( xp > -1 && xp < 1 ){
        //         return true;
        //     }else{
        //         return false;
        //     }
        // }
    public:
        double GetV0Instance(int index){
            return m_TrajectoryData->m_V0Instance[index];
        }
        double GetEta0Instance(int index){
            return m_TrajectoryData->m_Eta0Instance[index];
        }
        double GetBarrierState0Instance(int index){
            return m_TrajectoryData->m_BarrierState0Instance[index];
        }
        size_t GetNumberOfFrameX(int index){
            return m_TrajectoryData->m_TrajX[index].size();
        }
        double GetX(int index, size_t frameID){
            return m_TrajectoryData->m_TrajX[index][frameID];
        }
        std::vector<double>& GetX(int index){
            return m_TrajectoryData->m_TrajX[index];
        }
        size_t GetNumberOfFrameV(int index){
            return m_TrajectoryData->m_TrajV[index].size();
        }
        double GetV(int index, size_t frameID){
            return m_TrajectoryData->m_TrajV[index][frameID];
        }
        std::vector<double>& GetV(int index){
            return m_TrajectoryData->m_TrajV[index];
        }
        size_t GetNumberOfFrameF(int index){
            return m_TrajectoryData->m_TrajF[index].size();
        }
        double GetF(int index, size_t frameID){
            return m_TrajectoryData->m_TrajF[index][frameID];
        }
        std::vector<double>& GetF(int index){
            return m_TrajectoryData->m_TrajF[index];
        }
        size_t GetNumberOfFramePE(int index){
            return m_TrajectoryData->m_TrajPE[index].size();
        }
        double GetPE(int index, size_t frameID){
            return m_TrajectoryData->m_TrajPE[index][frameID];
        }
        std::vector<double>& GetPE(int index){
            return m_TrajectoryData->m_TrajPE[index];
        }
        size_t GetNumberOfFrameKE(int index){
            return m_TrajectoryData->m_TrajKE[index].size();
        }
        double GetKE(int index, size_t frameID){
            return m_TrajectoryData->m_TrajKE[index][frameID];
        }
        std::vector<double>& GetKE(int index){
            return m_TrajectoryData->m_TrajKE[index];
        }
        // eta
        size_t GetNumberOfFrameEta(int index){
            return m_TrajectoryData->m_TrajEta[index].size();
        }
        double GetEta(int index, size_t frameID){
            return m_TrajectoryData->m_TrajEta[index][frameID];
        }
        std::vector<double>& GetEta(int index){
            return m_TrajectoryData->m_TrajEta[index];
        }
        size_t GetNumberOfFrameState(int index){
            return m_TrajectoryData->m_TrajState[index].size();
        }
        double GetState(int index, size_t frameID){
            return m_TrajectoryData->m_TrajState[index][frameID];
        }
        std::vector<double>& GetState(int index){
            return m_TrajectoryData->m_TrajState[index];
        }
        size_t GetNumberOfFrameBarrier(int index){
            return m_TrajectoryData->m_TrajBarrier[index].size();
        }
        double GetBarrier(int index, size_t frameID){
            return m_TrajectoryData->m_TrajBarrier[index][frameID];
        }
        std::vector<double>& GetBarrier(int index){
            return m_TrajectoryData->m_TrajBarrier[index];
        }
        std::vector<double>& GetBarrierTick(int index){
            return m_TrajectoryData->m_BarrierTimeTick[index];
        }
    public:
        virtual bool Thin(int newFreq, int ncore){
            int thinFreq = newFreq / m_SaveFreq;
            if( thinFreq * m_SaveFreq != newFreq ){
                return false;
            }else{
                m_TrajectoryData->Thin(thinFreq, ncore);
                m_SaveFreq = newFreq;
                return true;
            }
        }
        size_t GetFirstPassageTimeCount(int index){
            return m_TrajectoryData->m_FirstPassageTime[index].size();
        }
        double GetFirstPassageTime(int index, size_t index2){
            return m_TrajectoryData->m_FirstPassageTime[index][index2];
        }
        double GetFirstPassageTimeAverage(int index){
            return STATISTICS::AVERAGE(&m_TrajectoryData->m_FirstPassageTime[index]);
        }        
        size_t GetFirstPassageTimeCount2(int index){
            return m_TrajectoryData->m_FirstPassageTime2[index].size();
        }
        double GetFirstPassageTime2(int index, size_t index2){
            return m_TrajectoryData->m_FirstPassageTime2[index][index2];
        }
        double GetFirstPassageTimeAverage2(int index){
            return STATISTICS::AVERAGE(&m_TrajectoryData->m_FirstPassageTime2[index]);
        }
        void GetInitialBarrierHeightHistogram(int index, HISTOGRAM::HISTOGRAM *histOut){
            HISTOGRAM::HISTOGRAM hist;
            hist.Allocate(histOut->min, histOut->max, histOut->size);
            hist.HistogramVector(&m_TrajectoryData->m_InitialBarrierList[index]);
            hist.NormalizePercent();
// hist.Print();
// printf("count %10zd %13.7f %d\n\n", m_TrajectoryData->m_InitialBarrierList[index].size(), m_TrajectoryData->m_CurrentState[index], index);
            for(int i=0;i<histOut->size;i++){
                histOut->count[i] += hist.percent[i]; // put percent into count, for convenience later
            }
        }
        void GetInstantBarrierHeightHistogram(int index, HISTOGRAM::HISTOGRAM *histOut){
            HISTOGRAM::HISTOGRAM hist;
            hist.Allocate(histOut->min, histOut->max, histOut->size);
            hist.HistogramVector(&m_TrajectoryData->m_InstantBarrierList[index]);
            hist.NormalizePercent();
            for(int i=0;i<histOut->size;i++){
                histOut->count[i] += hist.percent[i]; // put percent into count, for convenience later
            }
        }
        double GetTransitionPathTimeAverage(int index){
            return STATISTICS::AVERAGE(&m_TrajectoryData->m_TransitionPathTime[index]);
        }        
        double GetTransitionPathTimeAverage2(int index){
            return STATISTICS::AVERAGE(&m_TrajectoryData->m_TransitionPathTime2[index]);
        }
        std::vector<double>& GetTransitionPathTimeAll2(int index){
            return m_TrajectoryData->m_TransitionPathTime2[index];
        }
        double GetFirstPassageTimeCGAverage2(int index, int sizeIndex){
            return STATISTICS::AVERAGE(&m_TrajectoryData->m_FirstPassageTimeCGData[sizeIndex].m_FirstPassageTimeCG2[index]);
        }
    public:
        double GetStateCount(int index, int iState){
            return m_TrajectoryData->m_StateCount[iState][index];
        }
    }; // LangevinEquationSIMDMultiple

    class TrajectoryDataOverdamped : public TrajectoryData{
    public:
        __mSIMDd m_RandomForce;
    public:
        virtual void WriteRestartFile(StreamerWriter* writer){
            TrajectoryData::WriteRestartFile(writer);

            alignas(SIMDBYTESIZE) double rf[SIMDWIDTH];
            _mmSIMD_store_pd(rf, m_RandomForce);

            for(int i=0;i<SIMDWIDTH;i++){
                writer->WriteDouble(rf[i]);
            }
        };
        virtual void ReadRestartFile(StreamerReader* reader){
            TrajectoryData::ReadRestartFile(reader);
            
            alignas(SIMDBYTESIZE) double rf[SIMDWIDTH];
            for(int i=0;i<SIMDWIDTH;i++){
                rf[i] = reader->ReadDouble();
            }
            m_RandomForce = _mmSIMD_load_pd(rf);
        };
    };

    class OverdampedLangevinEquationSIMDMultiple : public LangevinEquationSIMDMultiple{
    protected:
        __mSIMDd m_TimeStepGammaMassInverser;
        TrajectoryDataOverdamped* m_TrajectoryDataOverdamped;
    protected:
        virtual void Allocate(){
            try{
                m_TrajectoryData = new TrajectoryDataOverdamped;
                m_TrajectoryDataOverdamped = static_cast<TrajectoryDataOverdamped*>(m_TrajectoryData);
            }catch ( std::bad_alloc& e){
                fprintf(stderr, "Error: failed to allocate trajectory data %s\n", e.what());
                exit(0);
            }
        }
    public:
        virtual void IntegrateInitializeVariables(){
            LangevinEquationSIMDMultiple::IntegrateInitializeVariables();

            double mGammaInverse = 1.0 / (m_M0 * m_Gamma0);
            double timeStepGammaMassInverser = m_TimeStep0 * mGammaInverse;
            m_TimeStepGammaMassInverser = _mmSIMD_set1_pd(timeStepGammaMassInverser);

            double varianceOfRandomForceIntegration = sqrt(m_kB0 * m_T0 * mGammaInverse * m_TimeStep0 * 0.5);
            m_VarianceOfRandomForceIntegration = _mmSIMD_set1_pd(varianceOfRandomForceIntegration);

            // never record V even if required by user ....
            m_SaveV = false;
        }
        virtual void IntegrateStep0(){
            LangevinEquationSIMDMultiple::IntegrateStep0();

            // R_n 
            m_TrajectoryDataOverdamped->m_RandomForce = _mmSIMD_mul_pd(
                m_VarianceOfRandomForceIntegration, 
                m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe());
        }
        virtual void IntegrateStep(size_t iStep){
            // R_n+1
            __mSIMDd rf = _mmSIMD_mul_pd(
                m_VarianceOfRandomForceIntegration, 
                m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe());

            // X_n+1 = X_n + dt^2/2 / m / gamma * F_n + dt/2 * sqrt(kb * T / m / gamma) * (R_n + R_n+1)
            m_TrajectoryDataOverdamped->m_X = _mmSIMD_fmadd_pd(m_TrajectoryDataOverdamped->m_F, m_TimeStepGammaMassInverser, m_TrajectoryDataOverdamped->m_X);
            // traj->m_X = _mm256_add_pd(traj->m_X, _mm256_mul_pd(traj->m_F, m_TimeStepGammaMassInverser));
            m_TrajectoryDataOverdamped->m_X = _mmSIMD_add_pd(m_TrajectoryDataOverdamped->m_X, rf);
            m_TrajectoryDataOverdamped->m_X = _mmSIMD_add_pd(m_TrajectoryDataOverdamped->m_X, m_TrajectoryDataOverdamped->m_RandomForce);

            // m_X = m_X + m_F * m_HalfTimeStepGammaMassInverser +
            //       rf + m_RandomForce; // (R_n + R_n+1)
            // R_n = R_n_1+1
            m_TrajectoryDataOverdamped->m_RandomForce = rf;

            UpdatePeriodicIndices(iStep);

            // even without transition, need to update current barrier when going to another barrier
            UpdateBarrierState(iStep);
            UpdateBarrierStateEtaParamter();

            // F_n+1
            //m_TrajectoryDataOverdamped->m_F = Force();
            EnergyAndForce();
        }
        virtual void IntegrateStepEta(size_t iStep){

            // R_n+1
            __mSIMDd rf = _mmSIMD_mul_pd(
                m_VarianceOfRandomForceIntegration, 
                m_GaussianRandomNumberSequenceSIMD.GetRandomNumberSafe());

            // X_n+1 = X_n + dt^2/2 / m / gamma * F_n + dt/2 * sqrt(kb * T / m / gamma) * (R_n + R_n+1)
            m_TrajectoryDataOverdamped->m_X = _mmSIMD_fmadd_pd(m_TrajectoryDataOverdamped->m_F, m_TimeStepGammaMassInverser, m_TrajectoryDataOverdamped->m_X);
            // traj->m_X = _mm256_add_pd(traj->m_X, _mm256_mul_pd(traj->m_F, m_TimeStepGammaMassInverser));
            m_TrajectoryDataOverdamped->m_X = _mmSIMD_add_pd(m_TrajectoryDataOverdamped->m_X, rf);
            m_TrajectoryDataOverdamped->m_X = _mmSIMD_add_pd(m_TrajectoryDataOverdamped->m_X, m_TrajectoryDataOverdamped->m_RandomForce);

            // m_X = m_X + m_F * m_HalfTimeStepGammaMassInverser +
            //       rf + m_RandomForce; // (R_n + R_n+1)
            // R_n = R_n_1+1
            m_TrajectoryDataOverdamped->m_RandomForce = rf;

            UpdatePeriodicIndices(iStep);

            // have to update the barrier here
            UpdateBarrierStateEta(iStep);
            UpdateBarrierStateEtaParamter();
            // F_n+1
            //m_TrajectoryDataOverdamped->m_F = Force();
            EnergyAndForce();
        }

    };

    class LangevinEquationRunnerMPI{
    public:
        text outputBase = "langevinEquation";

        double x0;
        double v0;
        double gamma;
        double kb;
        double temperature;
        double mass;
        double timestep;
        size_t numberOfStep;

        size_t statisticsFreq;
        size_t saveFreq;

        bool bRandomizeVelocity;
        bool bRandomizeEta;
        bool bRandomizeState;
        bool bSameInitialBarrier;

        double periodicLength;
        double barrierHeight[2];
        int initialState;
        double slopeHalfLength;
        double q;
        double tau;
        double eta0;

        bool bWaitingTimeDistribution;
        bool bFirstPassageTime;
        bool bTransitionPathTime;
        bool bFirstPassageTimeCG;
        std::vector<int> firstPassageTimeCGSize;

        bool bReadRestart;
        bool bWriteRestart;
        text sInputRestart;
        text sOutputRestart;

        int ntrial;
        size_t numberOfStepPrevious;

        bool writeEtaResult;
        double probHMin;
        double probHMax;
        double probHBin;
        std::vector<double> probTime;

        bool writeX;
        int samplePerOrder;
        int randomNumberSeedForSeed;

        bool bXCorr;
        bool bTimer;

        bool bGetV0Distribution;
        bool bGetEta0Distribution;
        bool bGetState0Distribution;

        bool bOverdamped;

        HISTOGRAM::HistogramSetting m_BarrierHistogramSetting;
        //
        MPIKERNEL::MPIKernel m_MPI;

        bool bEqualProbability;
        double equilibriumConstant;
    protected:
        int nrunPerNode;
        int nrun;
        int randomNumberSeedUsed;
        double numberOfTrajectorySet;
        double memoryUsage;
        double saveTimeStep;
    public:
        LangevinEquationRunnerMPI(MPI_Comm m_World){
            outputBase = "langevinEquation";
            x0 = 0.0;
            v0 = 0.0;
            gamma = 1.0;
            kb = 1.0;
            temperature = 1.0;
            mass = 1.0;
            timestep = 0.1;
            numberOfStep = 10000;

            statisticsFreq = 1;
            saveFreq = 1;

            bRandomizeVelocity = false;
            bRandomizeEta = false;
            bRandomizeState = true;

            periodicLength = 5.0;
            barrierHeight[0] = 2.0;
            barrierHeight[0] = 6.0;
            initialState = 0;
            slopeHalfLength = 1.0;
            q = 0.2;
            tau = 0.1;
            eta0 = 0.0;
            
            bWaitingTimeDistribution = false;
            bFirstPassageTime = false;
            bTransitionPathTime = false;
            bFirstPassageTimeCG = false;

            bReadRestart = false;
            bWriteRestart = false;
            sInputRestart = "inputRestart";
            sOutputRestart = "outputRestart";

            ntrial = 100000;
            numberOfStepPrevious = 0;

            writeEtaResult = false;
            probHMin = -100.0;
            probHMax = 100.0;
            probHBin = 0.1;
            probTime = {0.1, 1.0, 10.0, 100.0};

            writeX = false;
            samplePerOrder = 20;
            randomNumberSeedForSeed = -1;

            bXCorr = false;
            bTimer = false;

            bGetV0Distribution = false;
            bGetEta0Distribution = false;
            bGetState0Distribution = false;

            bOverdamped = false;

            m_BarrierHistogramSetting.Set(0.0, 10.0, 100);

            m_MPI.Initialize(m_World);

            nrunPerNode = 1;
            nrun = 1;
            randomNumberSeedUsed = -1;
            numberOfTrajectorySet = 0;
            memoryUsage = 0;
            saveTimeStep = 0.0;

            bEqualProbability = true;
            equilibriumConstant = 1.0;
        }
    public:
        void PrintParameter(FILE* fout = stderr){
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Overdamped           = %13s\n", text(bOverdamped).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Number of node       = %13d\n", m_MPI.GetNumberOfNode());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Time step            = %13.7f\n", timestep);
            fprintf(fout, "# Save frequency       = %13zd\n", saveFreq);
            fprintf(fout, "# Save time step       = %13.7f\n", saveTimeStep);
            fprintf(fout, "# Statistics freq      = %13zd\n", statisticsFreq);
            fprintf(fout, "# Number of step (pre) = %13zd\n", numberOfStepPrevious);
            fprintf(fout, "# Number of step       = %13zd\n", numberOfStep);
            fprintf(fout, "# Simulation time      = %13.6E\n", timestep * numberOfStep);
            fprintf(fout, "# Number of trial      = %13d\n", ntrial);
            fprintf(fout, "# Number of run/trial  = %13d\n", SIMDWIDTH);
            fprintf(fout, "# Number of run        = %13d\n", nrun);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Barrier height [0]   = %13.7f\n", barrierHeight[0]);
            fprintf(fout, "# Barrier height [1]   = %13.7f\n", barrierHeight[1]);
            fprintf(fout, "# Barrier periodicity  = %13.7f\n", periodicLength);
            fprintf(fout, "# slope half length    = %13.7f\n", slopeHalfLength);
            fprintf(fout, "# Equal probability    = %13s\n", text(bEqualProbability).c_str());
            fprintf(fout, "# Equilibrium constant = %13.7E\n", equilibriumConstant);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Mass                 = %13.7f\n", mass);
            fprintf(fout, "# Gamma                = %13.7f\n", gamma);
            fprintf(fout, "# Temperature          = %13.7f\n", temperature);
            fprintf(fout, "# kB                   = %13.7f\n", kb);
            fprintf(fout, "# Eta variance (Q)     = %13.7f\n", q);
            fprintf(fout, "# Eta tau              = %13.7f\n", tau);
            fprintf(fout, "# Sample per order     = %13d\n", samplePerOrder);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# X0                   = %13.7f\n", x0);
            if( bRandomizeVelocity ){
                fprintf(fout, "# V0                   = %13s\n", "Random");
            }else{
                fprintf(fout, "# V0                   = %13.7f\n", v0);
            }
            if( bRandomizeEta ){
                fprintf(fout, "# Eta0                 = %13s\n", "Random");
            }else{
                fprintf(fout, "# Eta0                 = %13.7f\n", eta0);
            }
            if( bRandomizeState ){
                fprintf(fout, "# State0               = %13s\n", "Random");
            }else{
                fprintf(fout, "# State0               = %13d\n", initialState);
            }
            fprintf(fout, "# Same initial barrier = %13s\n", text(bSameInitialBarrier).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Random number seed0  = %13d\n", randomNumberSeedForSeed);
            fprintf(fout, "# Random number seed   = %13d\n", randomNumberSeedUsed);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Number of traj(each) = %13.0f\n", numberOfTrajectorySet);
            fprintf(fout, "# Memory usage (each)  = %13s\n", text::FromByte2Readable(memoryUsage).c_str());
            fprintf(fout, "# Memory usage (total) = %13s\n", text::FromByte2Readable(memoryUsage * m_MPI.GetNumberOfNode()).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# MSD                  = %13s\n", "Mandatory");
            fprintf(fout, "# Waiting time         = %13s\n", text(bWaitingTimeDistribution).c_str());
            fprintf(fout, "# Write trajectories   = %13s\n", text(writeX).c_str());
            fprintf(fout, "# V0 Distribution      = %13s\n", text(bGetV0Distribution).c_str());
            fprintf(fout, "# Eta0 Distribution    = %13s\n", text(bGetEta0Distribution).c_str());
            fprintf(fout, "# State0 Distribution  = %13s\n", text(bGetState0Distribution).c_str());
            fprintf(fout, "# Eta result           = %13s\n", text(writeEtaResult).c_str());
            fprintf(fout, "# First passage time   = %13s\n", text(bFirstPassageTime).c_str());
            if( bFirstPassageTime ){
                fprintf(fout, "# Barrier hist hmin    = %13.7f\n", m_BarrierHistogramSetting.m_HistHMin);
                fprintf(fout, "# Barrier hist hmax    = %13.7f\n", m_BarrierHistogramSetting.m_HistHMax);
                fprintf(fout, "# Barrier hist nbin    = %13d\n", m_BarrierHistogramSetting.m_HistNBin);
            }
            fprintf(fout, "# Transtion path time  = %13s\n", text(bTransitionPathTime).c_str());
            fprintf(fout, "# First passage time CG= %13s\n", text(bFirstPassageTimeCG).c_str());
            if( bFirstPassageTimeCG ){
                for(int j=0;j<firstPassageTimeCGSize.size();j++){
                    fprintf(fout, "# FFT CG size %4d     = %13d\n", j, firstPassageTimeCGSize[j]);
                }
            }
            fprintf(fout, "# X correlation        = %13s\n", text(bXCorr).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Read restart         = %13s\n", text(bReadRestart).c_str());
            fprintf(fout, "# Write restart        = %13s\n", text(bWriteRestart).c_str());
            fprintf(fout, "# Input restart (base) = %13s\n", sInputRestart.c_str());
            fprintf(fout, "# Output restart (base)= %13s\n", sOutputRestart.c_str());
            fprintf(fout, "# ==============================================\n");
        }
        void Setup(
            text outputBase = "langevinEquation", 
            double x0 = 0.0, double v0 = 0.0, 
            double gamma = 1.0, double kb = 1.0, double temperature = 1.0, double mass = 1.0, 
            double timestep = 0.1, size_t numberOfStep = 10000, 
            size_t statisticsFreq = 100, size_t saveFreq = 1, 
            bool bRandomizeVelocity = false, bool bRandomizeEta = false, bool bRandomizeState = true, 
            bool bSameInitialBarrier = false,
            double periodicLength = 5.0, double barrierHeight0 = 1.0, double barrierHeight1 = 1.0, 
            int initialState = 0, double slopeHalfLength = 1.0, 
            double q = 0.2, double tau = 0.1, double eta0 = 0.0,
            bool bWaitingTimeDistribution = false, double bFirstPassageTime = false,
            bool bReadRestart = false, bool bWriteRestart = false, 
            text sInputRestart = "inputRestart", text sOutputRestart = "outputRestart", 
            int ntrial = 100000, size_t numberOfStepPrevious = 0, 
            bool writeEtaResult = false, double probHMin = -100.0, double probHMax = 100.0, 
            double probHBin = 0.1, std::vector<double> probTime = {0.1, 1.0, 10.0, 100.0},
            bool writeX = false, 
            int samplePerOrder = 20, int randomNumberSeedForSeed = -1,
            bool bXCorr = false, bool bTimer = false,
            bool bGetV0Distribution = false, bool bGetEta0Distribution = false, bool bGetState0Distribution = false,
            bool bOverdamped = false, 
            double barrierHistHMin = 0.0, double barrierHistHMax = 10.0, int barrierHistNBin = 100,
            bool bTransitionPathTime = false, bool bFirstPassageTimeCG = false, std::vector<int> firstPassageTimeCGSize = {},
            bool bEqualProbability = true, double equilibriumConstant = 1.0){
            // code
            this->outputBase = outputBase;
            this->x0 = x0;
            this->v0 = v0;
            this->gamma = gamma;
            this->kb = kb;
            this->temperature = temperature;
            this->mass = mass;
            this->timestep = timestep;
            this->numberOfStep = numberOfStep;

            this->statisticsFreq = statisticsFreq;
            this->saveFreq = saveFreq;

            this->bRandomizeVelocity = bRandomizeVelocity;
            this->bRandomizeEta = bRandomizeEta;
            this->bRandomizeState = bRandomizeState;
            this->bSameInitialBarrier = bSameInitialBarrier;

            this->periodicLength = periodicLength;
            this->barrierHeight[0] = barrierHeight0;
            this->barrierHeight[1] = barrierHeight1;
            this->initialState = initialState;
            this->slopeHalfLength = slopeHalfLength;
            this->q = q;
            this->tau = tau;
            this->eta0 = eta0;

            this->bWaitingTimeDistribution = bWaitingTimeDistribution;
            this->bFirstPassageTime = bFirstPassageTime;
            this->bTransitionPathTime = bTransitionPathTime;
            this->bFirstPassageTimeCG = bFirstPassageTimeCG;
            this->firstPassageTimeCGSize = firstPassageTimeCGSize;

            this->bReadRestart = bReadRestart;
            this->bWriteRestart = bWriteRestart;
            this->sInputRestart = sInputRestart;
            this->sOutputRestart = sOutputRestart;

            this->ntrial = ntrial;
            this->numberOfStepPrevious = numberOfStepPrevious;

            this->writeEtaResult = writeEtaResult;
            this->probHMin = probHMin;
            this->probHMax = probHMax;
            this->probHBin = probHBin;
            this->probTime = probTime;

            this->writeX = writeX;
            this->samplePerOrder = samplePerOrder;
            this->randomNumberSeedForSeed = randomNumberSeedForSeed;

            this->bXCorr = bXCorr;
            this->bTimer = bTimer;

            this->bGetV0Distribution = bGetV0Distribution;
            this->bGetEta0Distribution = bGetEta0Distribution;
            this->bGetState0Distribution = bGetState0Distribution;

            this->bOverdamped = bOverdamped;

            this->m_BarrierHistogramSetting.Set(barrierHistHMin, barrierHistHMax, barrierHistNBin);

            this->bEqualProbability = bEqualProbability;
            this->equilibriumConstant = equilibriumConstant;
        }
        void Run(){
            // code
            if( m_MPI.IsMasterNode() ) fprintf(stderr, "LangevinEquationMultipleRunner\n");

            size_t numberOfFrame = LangevinEquationSIMDMultiple::CalculateNumberOfFrame(numberOfStep + numberOfStepPrevious, saveFreq);
            if( numberOfFrame < 1000 ){
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "Error: numberOfFrame < 1000\n");
                return;
            }

            nrunPerNode = ntrial * SIMDWIDTH;
            nrun = nrunPerNode * m_MPI.GetNumberOfNode();

            std::vector<int> randomNumberSeedList = m_MPI.RandomSeed(randomNumberSeedForSeed);

            LangevinEquationSIMDMultiple* leOMP = 0;
            {
                if( bOverdamped ){
                    try{
                        leOMP = new OverdampedLangevinEquationSIMDMultiple();
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", m_MPI.GetNodeIndex(), e.what());
                        exit(0);
                    }
                }else{
                    try{
                        leOMP = new LangevinEquationSIMDMultiple();
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", m_MPI.GetNodeIndex(), e.what());
                        exit(0);
                    }
                }

                leOMP->m_X0 = x0;
                leOMP->m_V0 = v0;
                leOMP->m_Gamma0 = gamma;
                leOMP->m_kB0 = kb;
                leOMP->m_T0 = temperature;
                leOMP->m_M0 = mass;
                leOMP->m_TimeStep0 = timestep;
                leOMP->m_NumberOfStep = numberOfStep;

                leOMP->m_StatisticsFreq = statisticsFreq;
                leOMP->m_SaveFreq = saveFreq;

                leOMP->m_SaveX = true;
                if( writeX ){
                    leOMP->m_SaveV = true;
                    leOMP->m_SaveF = true;
                    leOMP->m_SaveE = true;
                    leOMP->m_SaveState = true;
                    leOMP->m_SaveBarrier = true;
                }
                if( writeEtaResult ){
                    leOMP->m_SaveEta = true;
                }

                leOMP->m_RandomizeVelocity = bRandomizeVelocity;
                leOMP->m_RandomizeEta = bRandomizeEta;
                leOMP->m_RandomizeState0 = bRandomizeState;
                leOMP->m_bSameInitialBarrier = bSameInitialBarrier;

                leOMP->m_PeriodicLength0 = periodicLength;
                leOMP->m_BarrierHeight0[0] = barrierHeight[0];
                leOMP->m_BarrierHeight0[1] = barrierHeight[1];
                leOMP->m_SlopeHalfLength0 = slopeHalfLength;
                leOMP->m_InitialState0 = initialState;
                leOMP->m_Q0 = q;
                leOMP->m_Tau0 = tau;
                leOMP->m_Eta0 = eta0;

                leOMP->m_bWaitingTimeDistribution = bWaitingTimeDistribution;
                leOMP->m_bFirstPassageTime = bFirstPassageTime;
                leOMP->m_bTransitionPathTime = bTransitionPathTime;
                leOMP->m_bFirstPassageTimeCG = bFirstPassageTimeCG;
                leOMP->m_FirstPassageTimeCGSize = firstPassageTimeCGSize;
                
                leOMP->m_ReadRestartFile = bReadRestart;
                leOMP->m_WriteRestartFile = bWriteRestart;
                leOMP->m_InputRestartFileName = sInputRestart;
                leOMP->m_OutputRestartFileName = sOutputRestart;

                leOMP->m_bEqualProbability = bEqualProbability;
                leOMP->m_EquilibriumConstant = equilibriumConstant;
                
                leOMP->m_Rank = m_MPI.GetNodeIndex();

                // this depends on tau 
                // randomNumberSeedUsed = leOMP->Initialize(randomNumberSeed);
                randomNumberSeedUsed = randomNumberSeedList[m_MPI.GetNodeIndex()];
                leOMP->Initialize(randomNumberSeedUsed);

                if( m_MPI.IsMasterNode() && ntrial != 1 && bReadRestart ){
                    fprintf(stderr, "Warning: Exact restart can only be produced of the first trial\n");
                    fprintf(stderr, "         Because the random number from the same seed is consumed differently\n");
                }
            }

            numberOfTrajectorySet = leOMP->CalculateNumberOfTrajectorySet();
            memoryUsage = leOMP->CalculateMemoryUsage();

            saveTimeStep = timestep * saveFreq;

            if( m_MPI.IsMasterNode() ){
                PrintParameter();
            }

            // int maxlag = maxLagTime / timestep;
            std::vector<size_t> interval = 
                MOLUTILITY::GenerateLog10Scale(
                    saveTimeStep, numberOfFrame, (size_t)1000, (size_t)samplePerOrder, true);
            if( interval.size() < 20 ){
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "Error: interval is less than 20 (%zd)\n", interval.size());
                SAFE_DELETE(leOMP);

                MPI_Finalize();
                exit(0);
            }else if( interval.size() < 100 ){
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "Warning: interval is less than 100 (%zd)\n", interval.size());
            }else{
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "Interval = %10zd\n", interval.size());
            }

            // msd (mandatory)
            std::vector<double> msd(interval.size());
            std::vector<double> msd0(interval.size());

            // state count
            double stateCount[2] = { 0, 0 };

            // xcorr
            std::vector<double> corr;
            std::vector<double> corr0;
            if( bXCorr ){
                corr.resize(interval.size());
                corr0.resize(interval.size());
            }

            // now the histogram is for eta
            HISTOGRAM::HistogramMultiple hist;
            std::vector<double> corrEta;
            std::vector<double> corrEta0;
            std::vector<size_t> probTimeStep;
            if( writeEtaResult ){
                corrEta.resize(interval.size());
                corrEta0.resize(interval.size());
                probTimeStep.resize(probTime.size());
                hist.Allocate(probTime.size(), probHMin, probHMax, probHBin);

                for(int i=0;i<probTime.size();i++){
                    probTimeStep[i] = probTime[i] / saveTimeStep;
                    if( probTimeStep[i] >= numberOfFrame ){
                        if( m_MPI.IsMasterNode() ) fprintf(stderr, "Error: eta prob time is too large\n");
                    }
                    probTime[i] = probTimeStep[i] * saveTimeStep;
                    if( m_MPI.IsMasterNode() ) fprintf(stderr, "Eta Prob at time %13.7f %10zd\n", 
                        probTime[i], probTimeStep[i]);
                }
            }


            // crossing event
            std::vector<double> waitingTimeList;

            // v0
            HISTOGRAM::HISTOGRAM histV0;
            if( bGetV0Distribution ){
                double velocitySigma = sqrt(kb * temperature / mass);
                histV0.Allocate(-velocitySigma * 5.0, velocitySigma * 5.0, 200);
            }

            // eta0
            HISTOGRAM::HISTOGRAM histEta0;
            if( bGetEta0Distribution ){
                histEta0.Allocate(-q * 5.0, q * 5.0, 200);
            }

            // s0
            HISTOGRAM::HISTOGRAM histState0;
            if( bGetState0Distribution ){
                histState0.Allocate(0.0, 2.0, 2);
            }

            // mean first passage time
            std::vector<double> firstPassageTime;
            std::vector<double> firstPassageTime2;
            HISTOGRAM::HISTOGRAM histInitialBarrierHeight;
            HISTOGRAM::HISTOGRAM histInstantBarrierHeight;
            if( bFirstPassageTime ){
                histInitialBarrierHeight.Allocate(m_BarrierHistogramSetting);
                histInstantBarrierHeight.Allocate(m_BarrierHistogramSetting);
            }

            // transitionPath time
            std::vector<double> transitionPathTime;
            std::vector<double> transitionPathTime2;

            // mean first passage time CG
            std::vector<std::vector<double> > firstPassageTimeCG2(firstPassageTimeCGSize.size());

            clsSystemTimeCPP timer;
            timer.AddEntry("total");

            clsSystemTimeCPP timer2;
            timer2.AddEntry("integrate");
            timer2.AddEntry("analysis");
            timer2.AddEntry("post");

            timer.Start();

            LoopTimer loopTimer;
            loopTimer.Start();

            for(int i=0;i<ntrial;i++){
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "##################################################################\n");
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "Trial %10d/%10d\n", i, ntrial);
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "##################################################################\n");
                
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "   Integrate\n");
                timer2.Start();
                bool result = leOMP->Integrate(i, bTimer && m_MPI.IsMasterNode());
                {
                    int iStatus = result ? 0 : 1;
                    int totalStatus = 0;
                    m_MPI.AllReduce(&iStatus, &totalStatus, 1, MPI_INT, MPI_SUM);
                    if( totalStatus != 0 ){
                        fprintf(stderr, "Error: Integrator failure\n");

                        SAFE_DELETE(leOMP);
                        exit(0);
                    }
                }
                timer2.StopToEntry("integrate");

                if( m_MPI.IsMasterNode() ) fprintf(stderr, "   Post Analysis\n");
                timer2.Start();
                for(int iChannel=0;iChannel<SIMDWIDTH;iChannel++){
                    // msd
                    MeanSquareDisplacementBruteForce(
                        &leOMP->GetX(iChannel), &msd0, &interval);
                    for(int t=0;t<interval.size();t++){
                        msd[t] += msd0[t];
                        if( isnan(msd0[t]) || isinf(msd0[t]) ){
                            fprintf(stderr, "Error: on node %10d, ill msd0\n", m_MPI.GetNodeIndex());
                        }
                    }
                    // x correlation
                    if( bXCorr ){ 
                        AutoCorrelationBruteForce(
                            &leOMP->GetX(iChannel), &corr0, &interval);
                        for(int t=0;t<interval.size();t++){
                            corr[t] += corr0[t];
                        }
                    }

                    // evaluate the PDF of eta
                    if( writeEtaResult ){
                        // histogram
                        for(int c=0;c<probTime.size();c++){
                            hist.AddValue(
                                leOMP->GetEta(iChannel, probTimeStep[c]), 1.0, c);
                        }
                        // correlation
                        std::vector<double> rawdata = leOMP->GetEta(iChannel);
                        STATISTICS::SUBTRACTAVERAGE(&rawdata);
                        AutoCorrelationBruteForce(
                            &rawdata, &corrEta0, &interval);

                        for(int t=0;t<interval.size();t++){
                            corrEta[t] += corrEta0[t];
                        }
                    }

                    // crossing event
                    // saved in time unit
                    if( bWaitingTimeDistribution ){
                        for(size_t t=1;t<leOMP->GetBarrierTick(iChannel).size();t++){
                            waitingTimeList.push_back(
                                (leOMP->GetBarrierTick(iChannel)[t] - 
                                leOMP->GetBarrierTick(iChannel)[t-1]) * timestep);
                        }
                    }

                    // V0
                    if( bGetV0Distribution ){
                        histV0.AddValue(leOMP->GetV0Instance(iChannel));
                    }

                    // Eta0
                    if( bGetEta0Distribution ){
                        histEta0.AddValue(leOMP->GetEta0Instance(iChannel));
                    }

                    // State0
                    if( bGetState0Distribution ){
                        histState0.AddValue(leOMP->GetBarrierState0Instance(iChannel));
                    }

                    if( writeX ){
                        char buffer[256];
                        sprintf(buffer, "%s_trajectory_%010d_%010d_%010d.log", 
                            outputBase.c_str(), i, iChannel, m_MPI.GetNodeIndex());
                        MOLUTILITY::findNextBackUpFileName(buffer);
                        FILE *fout = fopen(buffer, "w+");
                        PrintParameter(fout);
                        if( bOverdamped ){
                            for(size_t t=0;t<leOMP->GetNumberOfFrameX(iChannel);t++){
                                fprintf(fout, " %23.16f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n", 
                                    saveTimeStep * t, 
                                    leOMP->GetX(iChannel, t),
                                    0.0,
                                    leOMP->GetF(iChannel, t),
                                    leOMP->GetBarrier(iChannel, t),
                                    leOMP->GetPE(iChannel, t),
                                    0.0,
                                    leOMP->GetState(iChannel, t));
                            }
                        }else{
                            for(size_t t=0;t<leOMP->GetNumberOfFrameX(iChannel);t++){
                                fprintf(fout, " %23.16f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n", 
                                    saveTimeStep * t, 
                                    leOMP->GetX(iChannel, t),
                                    leOMP->GetV(iChannel, t),
                                    leOMP->GetF(iChannel, t),
                                    leOMP->GetBarrier(iChannel, t),
                                    leOMP->GetPE(iChannel, t),
                                    leOMP->GetKE(iChannel, t),
                                    leOMP->GetState(iChannel, t));
                            }
                        }
                        if( bWaitingTimeDistribution ){
                            fprintf(fout, "# event tick\n");
                            for(size_t t=0;t<leOMP->GetBarrierTick(iChannel).size();t++){
                                // internal ticks are in unit of timestep
                                fprintf(fout, "# %10zd %23.16E\n", t, leOMP->GetBarrierTick(iChannel)[t] * timestep);
                            }
                        }
                        fclose(fout);
                    }

                    if( bFirstPassageTime ){
                        firstPassageTime.push_back(leOMP->GetFirstPassageTimeAverage(iChannel));
                        firstPassageTime2.push_back(leOMP->GetFirstPassageTimeAverage2(iChannel));
                        leOMP->GetInitialBarrierHeightHistogram(iChannel, &histInitialBarrierHeight);
                        leOMP->GetInstantBarrierHeightHistogram(iChannel, &histInstantBarrierHeight);
                    }

                    if( bTransitionPathTime ){
                        transitionPathTime.push_back(leOMP->GetTransitionPathTimeAverage(iChannel));
                        transitionPathTime2.push_back(leOMP->GetTransitionPathTimeAverage2(iChannel));
                    }

                    if( bFirstPassageTimeCG ){
                        for(int jjj=0;jjj<firstPassageTimeCGSize.size();jjj++){
                            firstPassageTimeCG2[jjj].push_back(leOMP->GetFirstPassageTimeCGAverage2(iChannel, jjj));
                        }
                    }

                    double stateCount0[2] = { 0, 0 };
                    stateCount0[0] = leOMP->GetStateCount(iChannel, 0);
                    stateCount0[1] = leOMP->GetStateCount(iChannel, 1);

                    stateCount[0] += stateCount0[0];
                    stateCount[1] += stateCount0[1];
                }

                timer2.StopToEntry("analysis");

                loopTimer.Update(i, ntrial);
                if( m_MPI.IsMasterNode() ) fprintf(stderr, "   Done %10d/%10d %s\n", i + 1, ntrial, loopTimer.GetPChar());

                // synchronize
                m_MPI.Barrier();
            }

            m_MPI.Barrier();
            timer2.Start();

            if( bFirstPassageTime ){

                auto firstPassageTimeCombined = m_MPI.CombineVectorDoubleToMaster(&firstPassageTime);
                auto firstPassageTimeCombined2 = m_MPI.CombineVectorDoubleToMaster(&firstPassageTime2);
                m_MPI.SumVectorDouble(&histInitialBarrierHeight.count);
                m_MPI.SumVectorDouble(&histInstantBarrierHeight.count);

                if( m_MPI.IsMasterNode() ){
                    text filename = outputBase + "_firstPassageTime.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<firstPassageTimeCombined.size();index++){
                        writer.WriteDouble(firstPassageTimeCombined[index]);
                    }
                    writer.Close();
                }
                if( m_MPI.IsMasterNode() ){
                    text filename = outputBase + "_firstPassageTime2.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<firstPassageTimeCombined2.size();index++){
                        writer.WriteDouble(firstPassageTimeCombined2[index]);
                    }
                    writer.Close();
                }

                histInitialBarrierHeight.NormalizePercent();
                histInstantBarrierHeight.NormalizePercent();

                if( m_MPI.IsMasterNode() ){
                    {
                        text filename = outputBase + "_firstPassageTime.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&firstPassageTimeCombined));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&firstPassageTimeCombined));
                        fprintf(fout, "# count   %23zd\n", firstPassageTimeCombined.size());
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_firstPassageTime2.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&firstPassageTimeCombined2));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&firstPassageTimeCombined2));
                        fprintf(fout, "# count   %23zd\n", firstPassageTimeCombined2.size());
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_initialBarrierHeightHistogram.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        histInitialBarrierHeight.Print(fout);
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_instantBarrierHeightHistogram.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        histInstantBarrierHeight.Print(fout);
                        fclose(fout);
                    }
                }
            }

            {
                // state count
                double stateCountCombined[2] = { 0, 0 };
                m_MPI.AllReduce(stateCount, stateCountCombined, 2, MPI_DOUBLE, MPI_SUM);

                if( m_MPI.IsMasterNode() ){
                    text filename = outputBase + "_stateCount.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE* fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    fprintf(fout, "# state 0 %23.0f %13.7f\n", stateCountCombined[0], stateCountCombined[0] / (stateCountCombined[0] + stateCountCombined[1]));
                    fprintf(fout, "# state 1 %23.0f %13.7f\n", stateCountCombined[1], stateCountCombined[1] / (stateCountCombined[0] + stateCountCombined[1]));
                    fclose(fout);
                }
            }

            if( bTransitionPathTime ){

                auto transitionPathTimeCombined = m_MPI.CombineVectorDoubleToMaster(&transitionPathTime);
                auto transitionPathTimeCombined2 = m_MPI.CombineVectorDoubleToMaster(&transitionPathTime2);
                m_MPI.SumVectorDouble(&histInitialBarrierHeight.count);
                m_MPI.SumVectorDouble(&histInstantBarrierHeight.count);

                if( m_MPI.IsMasterNode() ){
                    text filename = outputBase + "_transitionPathTime.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<transitionPathTimeCombined.size();index++){
                        writer.WriteDouble(transitionPathTimeCombined[index]);
                    }
                    writer.Close();
                }
                if( m_MPI.IsMasterNode() ){
                    text filename = outputBase + "_transitionPathTime2.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<transitionPathTimeCombined2.size();index++){
                        writer.WriteDouble(transitionPathTimeCombined2[index]);
                    }
                    writer.Close();
                }

                if( m_MPI.IsMasterNode() ){
                    {
                        text filename = outputBase + "_transitionPathTime.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&transitionPathTimeCombined));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&transitionPathTimeCombined));
                        fprintf(fout, "# count   %23zd\n", transitionPathTimeCombined.size());
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_transitionPathTime2.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&transitionPathTimeCombined2));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&transitionPathTimeCombined2));
                        fprintf(fout, "# count   %23zd\n", transitionPathTimeCombined2.size());
                        fclose(fout);
                    }
                }
            }

            if( bFirstPassageTimeCG ){
                std::vector<std::vector<double > > firstPassageTimeCGCombined2(firstPassageTimeCGSize.size());
                for(int j=0;j<firstPassageTimeCGSize.size();j++){
                    firstPassageTimeCGCombined2[j] = m_MPI.CombineVectorDoubleToMaster(&firstPassageTimeCG2[j]);
                }

                if( m_MPI.IsMasterNode() ){
                    for(int j=0;j<firstPassageTimeCGSize.size();j++){
                        char buffer[256];
                        sprintf(buffer, "%s_firstPassageTimeCG2_%04d_%04d.dat", 
                            outputBase.c_str(), j, firstPassageTimeCGSize[j]);
                        text filename = buffer;
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        StreamerWriter writer;
                        writer.Open(filename, std::ios::binary);
                        for(int index=0;index<firstPassageTimeCGCombined2[j].size();index++){
                            writer.WriteDouble(firstPassageTimeCGCombined2[j][index]);
                        }
                        writer.Close();
                    }
                }

                if( m_MPI.IsMasterNode() ){
                    {
                        text filename = outputBase + "_firstPassageTimeCG2.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        for(int j=0;j<firstPassageTimeCGSize.size();j++){
                            fprintf(fout, "# %4d %4d\n", j, firstPassageTimeCGSize[j]);
                            fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&firstPassageTimeCGCombined2[j]));
                            fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&firstPassageTimeCGCombined2[j]));
                            fprintf(fout, "# count   %23zd\n", firstPassageTimeCGCombined2[j].size());
                        }
                        fclose(fout);
                    }
                }
            }

            // v0, eta0
            if( bGetV0Distribution ){
                // mpi sum
                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&histV0.count);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    histV0.NormalizePercent();

                    text filename = outputBase + "_V0Distribution.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    histV0.Print(fout);
                    fclose(fout);
                }
            }

            if( bGetEta0Distribution ){
                // mpi sum
                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&histEta0.count);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    histEta0.NormalizePercent();

                    text filename = outputBase + "_Eta0Distribution.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    histEta0.Print(fout);
                    fclose(fout);
                }
            }

            if( bGetState0Distribution ){
                // mpi sum
                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&histState0.count);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    histState0.NormalizePercent();

                    text filename = outputBase + "_State0Distribution.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    histState0.Print(fout);
                    fclose(fout);
                }
            }

            // now this is for eta
            if( writeEtaResult ){
                for(int i=0;i<probTime.size();i++){
                    m_MPI.Barrier();
                    m_MPI.SumVectorDouble(&hist.m_Data[i]);
                }

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    hist.NormalizePercent();

                    text filename = outputBase + "_EtaPDF.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int i=0;i<probTime.size();i++){
                        fprintf(fout, "# prob t = %13.7f %10zd\n", probTime[i], probTimeStep[i]);
                    }
                    hist.Print(fout);
                    fclose(fout);
                }

                // average corrEta
                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&corrEta);
                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    for(size_t t=0;t<interval.size();t++){
                        corrEta[t] /= double(nrun);
                    }

                    text filename = outputBase + "_EtaCorrelation.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int t=0;t<interval.size();t++){
                        fprintf(fout, "%13.7f %13.7f\n", saveTimeStep * interval[t], corrEta[t]);
                    }
                    fclose(fout);
                }
            }


            // x correlation
            if( bXCorr ) {
                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&corr);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    for(size_t t=0;t<interval.size();t++){
                        corr[t] /= double(nrun);
                    }

                    text filename = outputBase + "_correlation.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int t=0;t<interval.size();t++){
                        fprintf(fout, "%13.7f %13.7f\n", saveTimeStep * interval[t], corr[t]);
                    }
                    fclose(fout);
                }
            }

            // msd
            {
                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&msd);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    for(size_t t=0;t<interval.size();t++){
                        msd[t] /= double(nrun);
                    }

                    text filename = outputBase + "_msd.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int t=0;t<interval.size();t++){
                        fprintf(fout, "%14.7E %14.7E\n", saveTimeStep * interval[t], msd[t]);
                    }
                    fclose(fout);
                }
            }

            // crossing event
            if( bWaitingTimeDistribution ){
                double waitingTimeAvg = STATISTICS::SUM(&waitingTimeList);
                double waitingTimeAvg2 = STATISTICS::SUM2(&waitingTimeList);
                double waitingTimeMax = STATISTICS::MAX(&waitingTimeList);
                double waitingTimeListCount = waitingTimeList.size();

                {
                    double buffer[3], buffer2[3];
                    buffer[0] = waitingTimeAvg;
                    buffer[1] = waitingTimeAvg2;
                    buffer[2] = waitingTimeListCount;

                    m_MPI.Barrier();
                    m_MPI.AllReduce(buffer, buffer2, 3, MPI_DOUBLE, MPI_SUM);
                    waitingTimeAvg = buffer2[0];
                    waitingTimeAvg2 = buffer2[1];
                    waitingTimeListCount = buffer2[2];
                    waitingTimeAvg /= waitingTimeListCount;
                    waitingTimeAvg2 /= waitingTimeListCount;
                }
                
                {
                    double buffer, buffer2;
                    buffer = waitingTimeMax;
                    m_MPI.Barrier();
                    m_MPI.AllReduce(&buffer, &buffer2, 1, MPI_DOUBLE, MPI_MAX);
                    waitingTimeMax = buffer2;
                }

                double R = (waitingTimeAvg2 - waitingTimeAvg * waitingTimeAvg) / waitingTimeAvg / waitingTimeAvg;

                HISTOGRAM::HISTOGRAM histWaitingTime;
                int nbin = 100;
                if( waitingTimeListCount < 1000 ){
                    nbin = std::max<int>(10, (int)waitingTimeListCount / 100);
                }
                histWaitingTime.Allocate(0.0, waitingTimeMax + 0.1, nbin);
                histWaitingTime.HistogramVector(&waitingTimeList);

                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&histWaitingTime.count);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    histWaitingTime.NormalizePercent();

                    text filename = outputBase + "_waitingTime.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    fprintf(fout, "# R      = %13.7f\n", R);
                    fprintf(fout, "# Nwt    = %13zd\n", (size_t)waitingTimeListCount);
                    histWaitingTime.Print(fout);
                    fclose(fout);
                }
            }

            SAFE_DELETE(leOMP);
            timer2.StopToEntry("post");

            timer.StopToEntry("total");

            // time log
            {
                double timeTotal      = timer.GetTime("total");
                double timeIntegrate  = timer2.GetTime("integrate");
                double timeAnalysis   = timer2.GetTime("analysis");
                double timePost       = timer2.GetTime("post");
                double timeTotal2     = timeTotal * timeTotal;
                double timeIntegrate2 = timeIntegrate * timeIntegrate;
                double timeAnalysis2  = timeAnalysis * timeAnalysis;
                double timePost2      = timePost * timePost;


            
                std::vector<double> buffer(8);
                buffer[0] = timeTotal;
                buffer[1] = timeIntegrate;
                buffer[2] = timeAnalysis;
                buffer[3] = timePost;
                buffer[4] = timeTotal2;
                buffer[5] = timeIntegrate2;
                buffer[6] = timeAnalysis2;
                buffer[7] = timePost2;

                m_MPI.Barrier();
                m_MPI.SumVectorDouble(&buffer);

                m_MPI.Barrier();
                if( m_MPI.IsMasterNode() ){
                    timeTotal      = buffer[0] / double(m_MPI.GetNumberOfNode());
                    timeIntegrate  = buffer[1] / double(m_MPI.GetNumberOfNode());
                    timeAnalysis   = buffer[2] / double(m_MPI.GetNumberOfNode());
                    timePost       = buffer[3] / double(m_MPI.GetNumberOfNode());
                    timeTotal2     = buffer[4] / double(m_MPI.GetNumberOfNode());
                    timeIntegrate2 = buffer[5] / double(m_MPI.GetNumberOfNode());
                    timeAnalysis2  = buffer[6] / double(m_MPI.GetNumberOfNode());
                    timePost2      = buffer[7] / double(m_MPI.GetNumberOfNode());

                    timeTotal2 = sqrt(timeTotal2 - timeTotal * timeTotal);
                    timeIntegrate2 = sqrt(timeIntegrate2 - timeIntegrate * timeIntegrate);
                    timeAnalysis2 = sqrt(timeAnalysis2 - timeAnalysis * timeAnalysis);
                    timePost2 = sqrt(timePost2 - timePost * timePost);

                    text filename = outputBase + "_time.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE* fout = fopen(filename.c_str(), "w+");
                    for(int i=0;i<m_MPI.GetNumberOfNode();i++){
                        fprintf(fout, "Rank %10d Seed %10d\n", i, randomNumberSeedList[i]);
                    }
                    fprintf(fout, "Time consumed %s %s\n", 
                        timer2.Print(timeTotal, "%f %s").c_str(), timer2.Print(timeTotal2, "%f %s").c_str());
                    fprintf(fout, "Integrate     %s %s\n", 
                        timer2.Print(timeIntegrate, "%f %s").c_str(), timer2.Print(timeIntegrate2, "%f %s").c_str());
                    fprintf(fout, "Analysis      %s %s\n", 
                        timer2.Print(timeAnalysis, "%f %s").c_str(), timer2.Print(timeAnalysis2, "%f %s").c_str());
                    fprintf(fout, "Post          %s %s\n", 
                        timer2.Print(timePost, "%f %s").c_str(), timer2.Print(timePost2, "%f %s").c_str());
                    
                    fclose(fout);
                }
            }
            

        };
    };

    // ==========================================================================================
    // 
    // LangevinEquationRestartFile2MSD
    //
    // ==========================================================================================
    
    class LangevinEquationRestartFile2MSD{
    public:
        text outputBase = "langevinEquation";
        int ntrial;
        int samplePerOrder;
        bool bOverdamped;
        bool bReadRestart;
        bool bWriteRestart;
        text sInputRestart;
        text sOutputRestart;
        bool bWaitingTimeDistribution;
        bool bFirstPassageTime;
        bool bTransitionPathTime;
        bool bFirstPassageTimeCG;
        //
        int m_WorldSize;
    protected:
        int nrunPerNode;
        int nrun;
    public:
        LangevinEquationRestartFile2MSD(){
            outputBase = "langevinEquation";
            ntrial = 100000;

            bOverdamped = false;

            bReadRestart = false;
            bWriteRestart = false;
            sInputRestart = "inputRestart";
            sOutputRestart = "outputRestart";

            bWaitingTimeDistribution = false;
            bFirstPassageTime = false;
            bTransitionPathTime = false;
            bFirstPassageTimeCG = false;

            m_WorldSize = 1;
        }
    public:
        void Setup(
            text outputBase = "langevinEquation",
            int ntrial = 100000,
            int samplePerOrder = 20,
            bool bOverdamped = false, 
            bool bReadRestart = false, bool bWriteRestart = false, 
            text sInputRestart = "inputRestart", text sOutputRestart = "outputRestart", int worldSize = 1,
            bool bWaitingTimeDistribution = false, bool bFirstPassageTime = false, bool bTransitionPathTime = false, bool bFirstPassageTimeCG = false){
            // code
            this->outputBase = outputBase;
            this->ntrial = ntrial;
            this->bOverdamped = bOverdamped;
            this->bReadRestart = bReadRestart;
            this->bWriteRestart = bWriteRestart;
            this->sInputRestart = sInputRestart;
            this->sOutputRestart = sOutputRestart;
            this->m_WorldSize = worldSize;
            this->samplePerOrder = samplePerOrder;
            this->bWaitingTimeDistribution = bWaitingTimeDistribution;
            this->bFirstPassageTime = bFirstPassageTime;
            this->bTransitionPathTime = bTransitionPathTime;
            this->bFirstPassageTimeCG = bFirstPassageTimeCG;
        }
        void Convert(){            
            nrunPerNode = ntrial * SIMDWIDTH;
            nrun = nrunPerNode * m_WorldSize;

            std::vector<std::vector<double> > msd(nrun);

            std::vector<size_t> interval0;
            double saveTimeStep0;

            LoopTimer loopTimer;
            fprintf(stderr, "\n");
            loopTimer.Start();
            for(int rank=0;rank<m_WorldSize;rank++){
                fprintf(stderr, "\33[F%10d/%10d %s\n", rank, m_WorldSize, loopTimer.GetPChar());

                LangevinEquationSIMDMultiple* leOMP = 0;
                if( bOverdamped ){
                    try{
                        leOMP = new OverdampedLangevinEquationSIMDMultiple;
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", rank, e.what());
                        exit(0);
                    }
                }else{
                    try{
                        leOMP = new LangevinEquationSIMDMultiple;
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", rank, e.what());
                        exit(0);
                    }
                }
                leOMP->m_Rank = rank;
                leOMP->m_ReadRestartFile = bReadRestart;
                leOMP->m_WriteRestartFile = bWriteRestart;
                leOMP->m_InputRestartFileName = sInputRestart;
                leOMP->m_OutputRestartFileName = sOutputRestart;
                leOMP->m_bWaitingTimeDistribution = bWaitingTimeDistribution;
                leOMP->m_bFirstPassageTime = bFirstPassageTime;
                leOMP->m_bTransitionPathTime = bTransitionPathTime;
                leOMP->m_bFirstPassageTimeCG = bFirstPassageTimeCG;
                // this depends on tau 
                // randomNumberSeedUsed = leOMP->Initialize(randomNumberSeed);
                leOMP->Initialize(-1);

                for(int iTrial=0;iTrial<ntrial;iTrial++){
                    if( !leOMP->ReadRestartFileAll(iTrial) ){
                        fprintf(stderr, "Error: can't open restart file\n");
                    }

                    double saveTimeStep = leOMP->m_TimeStep0 * leOMP->m_SaveFreq;

                    size_t numberOfFrame = leOMP->CalculateNumberOfFrame(leOMP->m_NumberOfStep, leOMP->m_SaveFreq);

                    std::vector<size_t> interval = 
                        MOLUTILITY::GenerateLog10Scale(
                            saveTimeStep, numberOfFrame, (size_t)1000, (size_t)samplePerOrder, true);

                    interval0 = interval;
                    saveTimeStep0 = saveTimeStep;

                    for(int iChannel=0;iChannel<SIMDWIDTH;iChannel++){
                        int memoryID = rank * ntrial * SIMDWIDTH +
                                       iTrial * SIMDWIDTH +
                                       iChannel; 
                        msd[memoryID].resize(interval.size());

                        if( leOMP->GetX(iChannel).size() != numberOfFrame ){
                            fprintf(stderr, "Error: frame %10zd vs %10zd\n", 
                                leOMP->GetX(iChannel).size(), numberOfFrame);
                        }

                        MeanSquareDisplacementBruteForce(
                            &leOMP->GetX(iChannel), &msd[memoryID], &interval);
                    }
                }
                
                SAFE_DELETE(leOMP);

                loopTimer.Update(rank, m_WorldSize);
                fprintf(stderr, "\33[F%10d/%10d %s\n", rank + 1, m_WorldSize, loopTimer.GetPChar());
            }

            fprintf(stderr, "Calculate msd statistics\n");
            std::vector<double> msdAverage;
            std::vector<double> msdSD;
            if( msd.size() ){
                msdAverage.resize(msd[0].size());
                msdSD.resize(msd[0].size());
                for(int i=0;i<msd.size();i++){
                    for(int j=0;j<msdAverage.size();j++){
                        msdAverage[j] += msd[i][j];
                        msdSD[j] += msd[i][j] * msd[i][j];
                    }
                }
                for(int j=0;j<msdAverage.size();j++){
                    msdAverage[j] /= msd.size();
                    msdSD[j] /= msd.size();
                    msdSD[j] = sqrt(msdSD[j] - msdAverage[j] * msdAverage[j]);
                }
            }

            fprintf(stderr, "Write msd file\n");
            {
                text filename = outputBase + "_msd.log";
                MOLUTILITY::findNextBackUpFileName(filename.string());
                FILE* fout = fopen(filename.c_str(), "w+");
                for(int i=0;i<interval0.size();i++){
                    fprintf(fout, "%14.7E %14.7E %14.7E", 
                        saveTimeStep0 * interval0[i], 
                        msdAverage[i],
                        msdSD[i]);
                    for(int j=0;j<msd.size();j++){
                        fprintf(fout, " %13.7f", msd[j][i]);
                    }
                    fprintf(fout, "\n");
                }
                fclose(fout);
            }
        }
    };

    // ==========================================================================================
    // 
    // LangevinEquationRestartFile2Trajectory
    //
    // ==========================================================================================
    
    class LangevinEquationRestartFile2Trajectory{
    public:
        text outputBase = "langevinEquation";
        int ntrial;
        int samplePerOrder;
        bool bOverdamped;
        bool bReadRestart;
        bool bWriteRestart;
        text sInputRestart;
        text sOutputRestart;
        bool bWaitingTimeDistribution;
        bool bFirstPassageTime;
        bool bTransitionPathTime;
        bool bFirstPassageTimeCG;
        //
        int nodeID;
        int iTrial0;
        int iChannel0;
    protected:
        int nrunPerNode;
        int nrun;
    public:
        LangevinEquationRestartFile2Trajectory(){
            outputBase = "langevinEquation";
            ntrial = 100000;

            bOverdamped = false;

            bReadRestart = false;
            bWriteRestart = false;
            sInputRestart = "inputRestart";
            sOutputRestart = "outputRestart";

            bWaitingTimeDistribution = false;
            bFirstPassageTime = false;
            bTransitionPathTime = false;
            bFirstPassageTimeCG = false;

            nodeID = 0;
            iTrial0 = 0;
            iChannel0 = 0;
        }
    public:
        void Setup(
            text outputBase = "langevinEquation",
            int ntrial = 100000,
            int samplePerOrder = 20,
            bool bOverdamped = false, 
            bool bReadRestart = false, bool bWriteRestart = false, 
            text sInputRestart = "inputRestart", text sOutputRestart = "outputRestart",
            bool bWaitingTimeDistribution = false, double bFirstPassageTime = false,
            int nodeID = 0, int iTrial0 = 0, int iChannel0 = 0, bool bTransitionPathTime = false, bool bFirstPassageTimeCG = false){
            // code
            this->outputBase = outputBase;
            this->ntrial = ntrial;
            this->bOverdamped = bOverdamped;
            this->bReadRestart = bReadRestart;
            this->bWriteRestart = bWriteRestart;
            this->sInputRestart = sInputRestart;
            this->sOutputRestart = sOutputRestart;
            this->samplePerOrder = samplePerOrder;
            this->bWaitingTimeDistribution = bWaitingTimeDistribution;
            this->bFirstPassageTime = bFirstPassageTime;
            this->bTransitionPathTime = bTransitionPathTime;
            this->nodeID = nodeID;
            this->iTrial0 = iTrial0;
            this->iChannel0 = iChannel0;
            this->bFirstPassageTimeCG = bFirstPassageTimeCG;
        }
        void Convert(){            
            nrunPerNode = ntrial * SIMDWIDTH;

            double saveTimeStep0;

            {
                int rank = nodeID;

                LangevinEquationSIMDMultiple* leOMP = 0;
                if( bOverdamped ){
                    try{
                        leOMP = new OverdampedLangevinEquationSIMDMultiple;
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", rank, e.what());
                        exit(0);
                    }
                }else{
                    try{
                        leOMP = new LangevinEquationSIMDMultiple;
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", rank, e.what());
                        exit(0);
                    }
                }
                leOMP->m_Rank = rank;
                leOMP->m_ReadRestartFile = bReadRestart;
                leOMP->m_WriteRestartFile = bWriteRestart;
                leOMP->m_InputRestartFileName = sInputRestart;
                leOMP->m_OutputRestartFileName = sOutputRestart;
                leOMP->m_bWaitingTimeDistribution = bWaitingTimeDistribution;
                leOMP->m_bFirstPassageTime = bFirstPassageTime;
                leOMP->m_bTransitionPathTime = bTransitionPathTime;
                leOMP->m_bFirstPassageTimeCG = bFirstPassageTimeCG;
                // this depends on tau 
                // randomNumberSeedUsed = leOMP->Initialize(randomNumberSeed);
                leOMP->Initialize(-1);

                {
                    int iTrial = iTrial0;
                    if( !leOMP->ReadRestartFileAll(iTrial) ){
                        fprintf(stderr, "Error: can't open restart file\n");
                    }

                    double saveTimeStep = leOMP->m_TimeStep0 * leOMP->m_SaveFreq;

                    size_t numberOfFrame = leOMP->CalculateNumberOfFrame(leOMP->m_NumberOfStep, leOMP->m_SaveFreq);

                    {
                        char buffer[256];
                        sprintf(buffer, "%s_trajectory_%010d_%010d_%010d.log", 
                            outputBase.c_str(), iTrial0, iChannel0, nodeID);
                        MOLUTILITY::findNextBackUpFileName(buffer);
                        FILE *fout = fopen(buffer, "w+");
                        if( bOverdamped ){
                            for(size_t t=0;t<leOMP->GetNumberOfFrameX(iChannel0);t++){
                                fprintf(fout, " %23.16f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n", 
                                    saveTimeStep * t, 
                                    leOMP->GetX(iChannel0, t),
                                    0.0,
                                    leOMP->GetNumberOfFrameF(iChannel0) != 0 ? leOMP->GetF(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFrameBarrier(iChannel0) != 0 ? leOMP->GetBarrier(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFramePE(iChannel0) != 0 ? leOMP->GetPE(iChannel0, t) : 0.0,
                                    0.0,
                                    leOMP->GetNumberOfFrameState(iChannel0) != 0 ? leOMP->GetState(iChannel0, t) : 0.0);
                            }
                        }else{
                            for(size_t t=0;t<leOMP->GetNumberOfFrameX(iChannel0);t++){
                                fprintf(fout, " %23.16f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n", 
                                    saveTimeStep * t, 
                                    leOMP->GetX(iChannel0, t),
                                    leOMP->GetNumberOfFrameV(iChannel0) != 0 ? leOMP->GetV(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFrameF(iChannel0) != 0 ? leOMP->GetF(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFrameBarrier(iChannel0) != 0 ? leOMP->GetBarrier(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFramePE(iChannel0) != 0 ? leOMP->GetPE(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFrameKE(iChannel0) != 0 ? leOMP->GetKE(iChannel0, t) : 0.0,
                                    leOMP->GetNumberOfFrameState(iChannel0) != 0 ? leOMP->GetState(iChannel0, t) : 0.0);
                            }
                        }
                        if( bWaitingTimeDistribution ){
                            fprintf(fout, "# event tick\n");
                            for(size_t t=0;t<leOMP->GetBarrierTick(iChannel0).size();t++){
                                // internal ticks are in unit of timestep
                                fprintf(fout, "# %10zd %23.16E\n", t, leOMP->GetBarrierTick(iChannel0)[t] * leOMP->m_TimeStep0);
                            }
                        }
                        fclose(fout);
                    }
                }
                
                SAFE_DELETE(leOMP);

            }
        }
    };

    class LangevinEquationRunnerMPIPostProcess{
    public:
        text outputBase = "langevinEquation";

        double x0;
        double v0;
        double gamma;
        double kb;
        double temperature;
        double mass;
        double timestep;
        size_t numberOfStep;

        size_t statisticsFreq;
        size_t saveFreq;

        bool bRandomizeVelocity;
        bool bRandomizeEta;
        bool bRandomizeState;
        bool bSameInitialBarrier;

        double periodicLength;
        double barrierHeight[2];
        int initialState;
        double slopeHalfLength;
        double xc;
        double yc;

        double q;
        double tau;
        double eta0;

        bool bWaitingTimeDistribution;
        bool bFirstPassageTime;
        bool bTransitionPathTime;
        bool bFirstPassageTimeCG;
        std::vector<int> firstPassageTimeCGSize;

        bool bReadRestart;
        // bool bWriteRestart;
        text sInputRestart;
        // text sOutputRestart;

        int ntrial;
        size_t numberOfStepPrevious;

        bool writeEtaResult;
        double probHMin;
        double probHMax;
        double probHBin;
        std::vector<double> probTime;

        bool writeX;
        int samplePerOrder;
        int randomNumberSeedForSeed;

        bool bXCorr;
        bool bTimer;

        bool bGetV0Distribution;
        bool bGetEta0Distribution;
        bool bGetState0Distribution;

        bool bOverdamped;

        HISTOGRAM::HistogramSetting m_BarrierHistogramSetting;
        //
        // MPIKERNEL::MPIKernel m_MPI;
        int numberOfMPIThreads;

        bool bEqualProbability;
        double equilibriumConstant;
    protected:
        int nrunPerNode;
        int nrun;
        int randomNumberSeedUsed;
        double numberOfTrajectorySet;
        double memoryUsage;
        double saveTimeStep;
    public:
        LangevinEquationRunnerMPIPostProcess(){
            outputBase = "langevinEquation";
            x0 = 0.0;
            v0 = 0.0;
            gamma = 1.0;
            kb = 1.0;
            temperature = 1.0;
            mass = 1.0;
            timestep = 0.1;
            numberOfStep = 10000;

            statisticsFreq = 1;
            saveFreq = 1;

            bRandomizeVelocity = false;
            bRandomizeEta = false;
            bRandomizeState = true;

            periodicLength = 5.0;
            barrierHeight[0] = 2.0;
            barrierHeight[1] = 6.0;
            initialState = 0;
            slopeHalfLength = 1.0;
            xc = 2.5;
            yc = 2.0;

            q = 0.2;
            tau = 0.1;
            eta0 = 0.0;
            
            bWaitingTimeDistribution = false;
            bFirstPassageTime = false;
            bTransitionPathTime = false;
            bFirstPassageTimeCG = false;

            bReadRestart = false;
            // bWriteRestart = false;
            sInputRestart = "inputRestart";
            // sOutputRestart = "outputRestart";

            ntrial = 100000;
            numberOfStepPrevious = 0;

            writeEtaResult = false;
            probHMin = -100.0;
            probHMax = 100.0;
            probHBin = 0.1;
            probTime = {0.1, 1.0, 10.0, 100.0};

            writeX = false;
            samplePerOrder = 20;
            randomNumberSeedForSeed = -1;

            bXCorr = false;
            bTimer = false;

            bGetV0Distribution = false;
            bGetEta0Distribution = false;
            bGetState0Distribution = false;

            bOverdamped = false;

            m_BarrierHistogramSetting.Set(0.0, 10.0, 100);

            // m_MPI.Initialize(m_World);

            numberOfMPIThreads = 1;

            nrunPerNode = 1;
            nrun = 1;
            randomNumberSeedUsed = -1;
            numberOfTrajectorySet = 0;
            memoryUsage = 0;
            saveTimeStep = 0.0;

            bEqualProbability = true;
            equilibriumConstant = 1.0;
        }
    public:
        void PrintParameter(FILE* fout = stderr){
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Overdamped           = %13s\n", text(bOverdamped).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Number of node       = %13d\n", numberOfMPIThreads);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Time step            = %13.7f\n", timestep);
            fprintf(fout, "# Save frequency       = %13zd\n", saveFreq);
            fprintf(fout, "# Save time step       = %13.7f\n", saveTimeStep);
            fprintf(fout, "# Statistics freq      = %13zd\n", statisticsFreq);
            fprintf(fout, "# Number of step (pre) = %13zd\n", numberOfStepPrevious);
            fprintf(fout, "# Number of step       = %13zd\n", numberOfStep);
            fprintf(fout, "# Simulation time      = %13.6E\n", timestep * numberOfStep);
            fprintf(fout, "# Number of trial      = %13d\n", ntrial);
            fprintf(fout, "# Number of run/trial  = %13d\n", SIMDWIDTH);
            fprintf(fout, "# Number of run        = %13d\n", nrun);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Barrier height [0]   = %13.7f\n", barrierHeight[0]);
            fprintf(fout, "# Barrier height [1]   = %13.7f\n", barrierHeight[1]);
            fprintf(fout, "# Xc                   = %13.7f\n", xc);
            fprintf(fout, "# Yc                   = %13.7f\n", yc);
            fprintf(fout, "# Barrier periodicity  = %13.7f\n", periodicLength);
            fprintf(fout, "# Slope half length    = %13.7f\n", slopeHalfLength);
            fprintf(fout, "# Equal probability    = %13s\n", text(bEqualProbability).c_str());
            fprintf(fout, "# Equilibrium constant = %13.7E\n", equilibriumConstant);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Mass                 = %13.7f\n", mass);
            fprintf(fout, "# Gamma                = %13.7f\n", gamma);
            fprintf(fout, "# Temperature          = %13.7f\n", temperature);
            fprintf(fout, "# kB                   = %13.7f\n", kb);
            fprintf(fout, "# Eta variance (Q)     = %13.7f\n", q);
            fprintf(fout, "# Eta tau              = %13.7f\n", tau);
            fprintf(fout, "# Sample per order     = %13d\n", samplePerOrder);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# X0                   = %13.7f\n", x0);
            if( bRandomizeVelocity ){
                fprintf(fout, "# V0                   = %13s\n", "Random");
            }else{
                fprintf(fout, "# V0                   = %13.7f\n", v0);
            }
            if( bRandomizeEta ){
                fprintf(fout, "# Eta0                 = %13s\n", "Random");
            }else{
                fprintf(fout, "# Eta0                 = %13.7f\n", eta0);
            }
            if( bRandomizeState ){
                fprintf(fout, "# State0               = %13s\n", "Random");
            }else{
                fprintf(fout, "# State0               = %13d\n", initialState);
            }
            fprintf(fout, "# Same initial barrier = %13s\n", text(bSameInitialBarrier).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Random number seed0  = %13d\n", randomNumberSeedForSeed);
            fprintf(fout, "# Random number seed   = %13d\n", randomNumberSeedUsed);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Number of traj(each) = %13.0f\n", numberOfTrajectorySet);
            fprintf(fout, "# Memory usage (each)  = %13s\n", text::FromByte2Readable(memoryUsage).c_str());
            fprintf(fout, "# Memory usage (total) = %13s\n", text::FromByte2Readable(memoryUsage * numberOfMPIThreads).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# MSD                  = %13s\n", "Mandatory");
            fprintf(fout, "# Waiting time         = %13s\n", text(bWaitingTimeDistribution).c_str());
            fprintf(fout, "# Write trajectories   = %13s\n", text(writeX).c_str());
            fprintf(fout, "# V0 Distribution      = %13s\n", text(bGetV0Distribution).c_str());
            fprintf(fout, "# Eta0 Distribution    = %13s\n", text(bGetEta0Distribution).c_str());
            fprintf(fout, "# State0 Distribution  = %13s\n", text(bGetState0Distribution).c_str());
            fprintf(fout, "# Eta result           = %13s\n", text(writeEtaResult).c_str());
            fprintf(fout, "# First passage time   = %13s\n", text(bFirstPassageTime).c_str());
            if( bFirstPassageTime ){
                fprintf(fout, "# Barrier hist hmin    = %13.7f\n", m_BarrierHistogramSetting.m_HistHMin);
                fprintf(fout, "# Barrier hist hmax    = %13.7f\n", m_BarrierHistogramSetting.m_HistHMax);
                fprintf(fout, "# Barrier hist nbin    = %13d\n", m_BarrierHistogramSetting.m_HistNBin);
            }
            fprintf(fout, "# Transtion path time  = %13s\n", text(bTransitionPathTime).c_str());
            fprintf(fout, "# First passage time CG= %13s\n", text(bFirstPassageTimeCG).c_str());
            if( bFirstPassageTimeCG ){
                for(int j=0;j<firstPassageTimeCGSize.size();j++){
                    fprintf(fout, "# FFT CG size %4d     = %13d\n", j, firstPassageTimeCGSize[j]);
                }
            }
            fprintf(fout, "# X correlation        = %13s\n", text(bXCorr).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Read restart         = %13s\n", text(bReadRestart).c_str());
            fprintf(fout, "# Write restart        = %13s\n", "false");
            fprintf(fout, "# Input restart (base) = %13s\n", sInputRestart.c_str());
            fprintf(fout, "# Output restart (base)= %13s\n", "none");
            fprintf(fout, "# ==============================================\n");
        }
        void Setup(
            text outputBase = "langevinEquation", 
            double x0 = 0.0, double v0 = 0.0, 
            double gamma = 1.0, double kb = 1.0, double temperature = 1.0, double mass = 1.0, 
            double timestep = 0.1, size_t numberOfStep = 10000, 
            size_t statisticsFreq = 100, size_t saveFreq = 1, 
            bool bRandomizeVelocity = false, bool bRandomizeEta = false, bool bRandomizeState = true, 
            bool bSameInitialBarrier = false,
            double periodicLength = 5.0, 
            double barrierHeightBL = 1.0, double barrierHeightBH = 1.0, 
            int initialState = 0, 
            double slopeHalfLengthB = 1.0, 
            double q = 0.2, double tau = 0.1, double eta0 = 0.0,
            bool bWaitingTimeDistribution = false, double bFirstPassageTime = false,
            bool bReadRestart = false, // bool bWriteRestart = false, 
            text sInputRestart = "inputRestart", // text sOutputRestart = "outputRestart", 
            int ntrial = 100000, size_t numberOfStepPrevious = 0, 
            bool writeEtaResult = false, double probHMin = -100.0, double probHMax = 100.0, 
            double probHBin = 0.1, std::vector<double> probTime = {0.1, 1.0, 10.0, 100.0},
            bool writeX = false, 
            int samplePerOrder = 20, int randomNumberSeedForSeed = -1,
            bool bXCorr = false, bool bTimer = false,
            bool bGetV0Distribution = false, bool bGetEta0Distribution = false, bool bGetState0Distribution = false,
            bool bOverdamped = false, 
            double barrierHistHMin = 0.0, double barrierHistHMax = 10.0, int barrierHistNBin = 100,
            bool bTransitionPathTime = false, bool bFirstPassageTimeCG = false, std::vector<int> firstPassageTimeCGSize = {},
            bool bEqualProbability = true, double equilibriumConstant = 1.0, int numberOfMPIThreads = 1){
            // code
            this->outputBase = outputBase;
            this->x0 = x0;
            this->v0 = v0;
            this->gamma = gamma;
            this->kb = kb;
            this->temperature = temperature;
            this->mass = mass;
            this->timestep = timestep;
            this->numberOfStep = numberOfStep;

            this->statisticsFreq = statisticsFreq;
            this->saveFreq = saveFreq;

            this->bRandomizeVelocity = bRandomizeVelocity;
            this->bRandomizeEta = bRandomizeEta;
            this->bRandomizeState = bRandomizeState;
            this->bSameInitialBarrier = bSameInitialBarrier;

            this->periodicLength = periodicLength;
            this->barrierHeight[0] = barrierHeightBL;
            this->barrierHeight[1] = barrierHeightBH;
            this->initialState = initialState;
            this->slopeHalfLength = slopeHalfLengthB;
            
            this->q = q;
            this->tau = tau;
            this->eta0 = eta0;

            this->bWaitingTimeDistribution = bWaitingTimeDistribution;
            this->bFirstPassageTime = bFirstPassageTime;
            this->bTransitionPathTime = bTransitionPathTime;
            this->bFirstPassageTimeCG = bFirstPassageTimeCG;
            this->firstPassageTimeCGSize = firstPassageTimeCGSize;

            this->bReadRestart = bReadRestart;
            // this->bWriteRestart = bWriteRestart;
            this->sInputRestart = sInputRestart;
            // this->sOutputRestart = sOutputRestart;

            this->ntrial = ntrial;
            this->numberOfStepPrevious = numberOfStepPrevious;

            this->writeEtaResult = writeEtaResult;
            this->probHMin = probHMin;
            this->probHMax = probHMax;
            this->probHBin = probHBin;
            this->probTime = probTime;

            this->writeX = writeX;
            this->samplePerOrder = samplePerOrder;
            this->randomNumberSeedForSeed = randomNumberSeedForSeed;

            this->bXCorr = bXCorr;
            this->bTimer = bTimer;

            this->bGetV0Distribution = bGetV0Distribution;
            this->bGetEta0Distribution = bGetEta0Distribution;
            this->bGetState0Distribution = bGetState0Distribution;

            this->bOverdamped = bOverdamped;

            this->m_BarrierHistogramSetting.Set(barrierHistHMin, barrierHistHMax, barrierHistNBin);

            this->bEqualProbability = bEqualProbability;
            this->equilibriumConstant = equilibriumConstant;

            this->numberOfMPIThreads = numberOfMPIThreads;
        }
        void Run(){
            // code
            fprintf(stderr, "LangevinEquationMultipleRunner\n");

            size_t numberOfFrame = LangevinEquationSIMDMultiple::CalculateNumberOfFrame(numberOfStep + numberOfStepPrevious, saveFreq);
            if( numberOfFrame < 1000 ){
                fprintf(stderr, "Error: numberOfFrame < 1000\n");
                return;
            }

            nrunPerNode = ntrial * SIMDWIDTH;
            nrun = nrunPerNode * numberOfMPIThreads;

            // std::vector<int> randomNumberSeedList = m_MPI.RandomSeed(randomNumberSeedForSeed);

            LangevinEquationSIMDMultiple* leOMP = 0;
            {
                if( bOverdamped ){
                    try{
                        leOMP = new OverdampedLangevinEquationSIMDMultiple();
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", 0, e.what());
                        exit(0);
                    }
                }else{
                    try{
                        leOMP = new LangevinEquationSIMDMultiple();
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", 0, e.what());
                        exit(0);
                    }
                }

                leOMP->m_X0 = x0;
                leOMP->m_V0 = v0;
                leOMP->m_Gamma0 = gamma;
                leOMP->m_kB0 = kb;
                leOMP->m_T0 = temperature;
                leOMP->m_M0 = mass;
                leOMP->m_TimeStep0 = timestep;
                leOMP->m_NumberOfStep = numberOfStep;

                leOMP->m_StatisticsFreq = statisticsFreq;
                leOMP->m_SaveFreq = saveFreq;

                leOMP->m_SaveX = true;
                if( writeX ){
                    leOMP->m_SaveV = true;
                    leOMP->m_SaveF = true;
                    leOMP->m_SaveE = true;
                    leOMP->m_SaveState = true;
                    leOMP->m_SaveBarrier = true;
                }
                if( writeEtaResult ){
                    leOMP->m_SaveEta = true;
                }

                leOMP->m_RandomizeVelocity = bRandomizeVelocity;
                leOMP->m_RandomizeEta = bRandomizeEta;
                leOMP->m_RandomizeState0 = bRandomizeState;
                leOMP->m_bSameInitialBarrier = bSameInitialBarrier;

                leOMP->m_PeriodicLength0 = periodicLength;
                leOMP->m_BarrierHeight0[0] = barrierHeight[0];
                leOMP->m_BarrierHeight0[1] = barrierHeight[1];
                leOMP->m_SlopeHalfLength0 = slopeHalfLength;

                leOMP->m_InitialState0 = initialState;
                leOMP->m_Q0 = q;
                leOMP->m_Tau0 = tau;
                leOMP->m_Eta0 = eta0;

                leOMP->m_bWaitingTimeDistribution = bWaitingTimeDistribution;
                leOMP->m_bFirstPassageTime = bFirstPassageTime;
                leOMP->m_bTransitionPathTime = bTransitionPathTime;
                leOMP->m_bFirstPassageTimeCG = bFirstPassageTimeCG;
                leOMP->m_FirstPassageTimeCGSize = firstPassageTimeCGSize;
                
                leOMP->m_ReadRestartFile = bReadRestart;
                leOMP->m_WriteRestartFile = false;
                leOMP->m_InputRestartFileName = sInputRestart;
                leOMP->m_OutputRestartFileName = "none";

                leOMP->m_bEqualProbability = bEqualProbability;
                leOMP->m_EquilibriumConstant = equilibriumConstant;
                
                leOMP->m_Rank = 0;

                // this depends on tau 
                // randomNumberSeedUsed = leOMP->Initialize(randomNumberSeed);
                randomNumberSeedUsed = -1;
                leOMP->Initialize(randomNumberSeedUsed);

                // if( m_MPI.IsMasterNode() && ntrial != 1 && bReadRestart ){
                //     fprintf(stderr, "Warning: Exact restart can only be produced of the first trial\n");
                //     fprintf(stderr, "         Because the random number from the same seed is consumed differently\n");
                // }
            }

            numberOfTrajectorySet = leOMP->CalculateNumberOfTrajectorySet();
            memoryUsage = leOMP->CalculateMemoryUsage();

            saveTimeStep = timestep * saveFreq;

            // if( m_MPI.IsMasterNode() ){
                PrintParameter();
            // }

            // int maxlag = maxLagTime / timestep;
            std::vector<size_t> interval = 
                MOLUTILITY::GenerateLog10Scale(
                    saveTimeStep, numberOfFrame, (size_t)1000, (size_t)samplePerOrder, true);
            if( interval.size() < 20 ){
                fprintf(stderr, "Error: interval is less than 20 (%zd)\n", interval.size());
                SAFE_DELETE(leOMP);

                MPI_Finalize();
                exit(0);
            }else if( interval.size() < 100 ){
                fprintf(stderr, "Warning: interval is less than 100 (%zd)\n", interval.size());
            }else{
                fprintf(stderr, "Interval = %10zd\n", interval.size());
            }

            // msd (mandatory)
            std::vector<std::vector<double> > msd;
            std::vector<std::vector<double> > msd0;
            VectorAuxiliary::AllocateVector2D(&msd, numberOfMPIThreads, interval.size());
            VectorAuxiliary::AllocateVector2D(&msd0, numberOfMPIThreads, interval.size());

            // state count
            std::vector<double> stateCount[2];
            stateCount[0].resize(numberOfMPIThreads);
            stateCount[1].resize(numberOfMPIThreads);

            // xcorr
            std::vector<std::vector<double> > corr;
            std::vector<std::vector<double> > corr0;
            if( bXCorr ){
                VectorAuxiliary::AllocateVector2D(&corr, numberOfMPIThreads, interval.size());
                VectorAuxiliary::AllocateVector2D(&corr0, numberOfMPIThreads, interval.size());
            }

            // now the histogram is for eta
            std::vector<HISTOGRAM::HistogramMultiple> hist;
            std::vector<std::vector<double> > corrEta;
            std::vector<std::vector<double> > corrEta0;
            std::vector<size_t> probTimeStep(interval.size());
            if( writeEtaResult ){
                VectorAuxiliary::AllocateVector2D(&corrEta, numberOfMPIThreads, interval.size());
                VectorAuxiliary::AllocateVector2D(&corrEta0, numberOfMPIThreads, interval.size());

                hist.resize(numberOfMPIThreads);

                for(int j=0;j<numberOfMPIThreads;j++){
                    hist[j].Allocate(probTime.size(), probHMin, probHMax, probHBin);
                }

                for(int i=0;i<probTime.size();i++){
                    probTimeStep[i] = probTime[i] / saveTimeStep;
                    if( probTimeStep[i] >= numberOfFrame ){
                        fprintf(stderr, "Error: eta prob time is too large\n");
                    }
                    probTime[i] = probTimeStep[i] * saveTimeStep;
                    fprintf(stderr, "Eta Prob at time %13.7f %10zd\n", 
                        probTime[i], probTimeStep[i]);
                }
            }


            // crossing event
            std::vector<std::vector<double> > waitingTimeList(numberOfMPIThreads);

            // v0
            std::vector<HISTOGRAM::HISTOGRAM> histV0;
            if( bGetV0Distribution ){
                double velocitySigma = sqrt(kb * temperature / mass);
                histV0.resize(numberOfMPIThreads);
                for(int i=0;i<numberOfMPIThreads;i++){
                    histV0[i].Allocate(-velocitySigma * 5.0, velocitySigma * 5.0, 200);
                }
            }

            // eta0
            std::vector<HISTOGRAM::HISTOGRAM> histEta0;
            if( bGetEta0Distribution ){
                histEta0.resize(numberOfMPIThreads);
                for(int i=0;i<numberOfMPIThreads;i++){
                    histEta0[i].Allocate(-q * 5.0, q * 5.0, 200);
                }
            }

            // s0
            std::vector<HISTOGRAM::HISTOGRAM> histState0;
            if( bGetState0Distribution ){
                histState0.resize(numberOfMPIThreads);
                for(int i=0;i<numberOfMPIThreads;i++){
                    histState0[i].Allocate(0.0, 2.0, 2);
                }
            }

            // mean first passage time
            std::vector<std::vector<double> > firstPassageTime;
            std::vector<std::vector<double> > firstPassageTime2;
            std::vector<HISTOGRAM::HISTOGRAM> histInitialBarrierHeight;
            std::vector<HISTOGRAM::HISTOGRAM> histInstantBarrierHeight;
            if( bFirstPassageTime ){
                firstPassageTime.resize(numberOfMPIThreads);
                firstPassageTime2.resize(numberOfMPIThreads);
                histInitialBarrierHeight.resize(numberOfMPIThreads);
                histInstantBarrierHeight.resize(numberOfMPIThreads);
                for(int i=0;i<numberOfMPIThreads;i++){
                    histInitialBarrierHeight[i].Allocate(m_BarrierHistogramSetting);
                    histInstantBarrierHeight[i].Allocate(m_BarrierHistogramSetting);
                }
            }

            // transitionPath time
            std::vector<std::vector<double> > transitionPathTime(numberOfMPIThreads);
            std::vector<std::vector<double> > transitionPathTime2(numberOfMPIThreads);
            std::vector<double> transitionPathTime2All;

            // mean first passage time CG
            std::vector<std::vector<std::vector<double> > > firstPassageTimeCG2;
            VectorAuxiliary::AllocateVector2D(&firstPassageTimeCG2, numberOfMPIThreads, firstPassageTimeCGSize.size());

            clsSystemTimeCPP timer;
            timer.AddEntry("total");

            clsSystemTimeCPP timer2;
            timer2.AddEntry("integrate");
            timer2.AddEntry("analysis");
            timer2.AddEntry("post");

            timer.Start();

            LoopTimer loopTimer;
            loopTimer.Start();

            for(int i=0;i<ntrial;i++){
                fprintf(stderr, "##################################################################\n");
                fprintf(stderr, "Trial %10d/%10d\n", i, ntrial);
                fprintf(stderr, "##################################################################\n");
                
                for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                    leOMP->m_Rank = iMPI;
                    fprintf(stderr, "   Read          %10d\n", iMPI);
                    timer2.Start();
                    if( !leOMP->ReadRestartFileAll(i, false) ){
                        fprintf(stderr, "Error: can't open restart file\n");
                        SAFE_DELETE(leOMP);
                        exit(0);
                    }
                    timer2.StopToEntry("integrate");

                    fprintf(stderr, "   Post Analysis %10d\n", iMPI);
                    timer2.Start();

                    for(int iChannel=0;iChannel<SIMDWIDTH;iChannel++){
                        // msd
                        MeanSquareDisplacementBruteForce(
                            &leOMP->GetX(iChannel), &msd0[iMPI], &interval);

                        for(int t=0;t<interval.size();t++){
                            msd[iMPI][t] += msd0[iMPI][t];
                            if( isnan(msd0[iMPI][t]) || isinf(msd0[iMPI][t]) ){
                                fprintf(stderr, "Error: on node %10d, ill msd0\n", iMPI);
                            }
                        }
                        // x correlation
                        if( bXCorr ){ 
                            AutoCorrelationBruteForce(
                                &leOMP->GetX(iChannel), &corr0[iMPI], &interval);
                            for(int t=0;t<interval.size();t++){
                                corr[iMPI][t] += corr0[iMPI][t];
                            }
                        }

                        // evaluate the PDF of eta
                        if( writeEtaResult ){
                            // histogram
                            for(int c=0;c<probTime.size();c++){
                                hist[iMPI].AddValue(
                                    leOMP->GetEta(iChannel, probTimeStep[c]), 1.0, c);
                            }
                            // correlation
                            std::vector<double> rawdata = leOMP->GetEta(iChannel);
                            STATISTICS::SUBTRACTAVERAGE(&rawdata);
                            AutoCorrelationBruteForce(
                                &rawdata, &corrEta0[iMPI], &interval);

                            for(int t=0;t<interval.size();t++){
                                corrEta[iMPI][t] += corrEta0[iMPI][t];
                            }
                        }

                        // crossing event
                        // saved in time unit
                        if( bWaitingTimeDistribution ){
                            for(size_t t=1;t<leOMP->GetBarrierTick(iChannel).size();t++){
                                waitingTimeList[iMPI].push_back(
                                    (leOMP->GetBarrierTick(iChannel)[t] - 
                                    leOMP->GetBarrierTick(iChannel)[t-1]) * timestep);
                            }
                        }

                        // V0
                        if( bGetV0Distribution ){
                            histV0[iMPI].AddValue(leOMP->GetV0Instance(iChannel));
                        }

                        // Eta0
                        if( bGetEta0Distribution ){
                            histEta0[iMPI].AddValue(leOMP->GetEta0Instance(iChannel));
                        }

                        // State0
                        if( bGetState0Distribution ){
                            histState0[iMPI].AddValue(leOMP->GetBarrierState0Instance(iChannel));
                        }

                        if( writeX ){
                            char buffer[256];
                            sprintf(buffer, "%s_trajectory_%010d_%010d_%010d.log", 
                                outputBase.c_str(), i, iChannel, iMPI);
                            MOLUTILITY::findNextBackUpFileName(buffer);
                            FILE *fout = fopen(buffer, "w+");
                            PrintParameter(fout);
                            if( bOverdamped ){
                                for(size_t t=0;t<leOMP->GetNumberOfFrameX(iChannel);t++){
                                    fprintf(fout, " %23.16f %23.16f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n", 
                                        saveTimeStep * t, 
                                        leOMP->GetX(iChannel, t),
                                        0.0,
                                        leOMP->GetNumberOfFrameF(iChannel) != 0 ? leOMP->GetF(iChannel, t) : 0.0,
                                        leOMP->GetNumberOfFrameBarrier(iChannel) != 0 ? leOMP->GetBarrier(iChannel, t) : 0.0,
                                        leOMP->GetNumberOfFramePE(iChannel) != 0 ? leOMP->GetPE(iChannel, t) : 0.0,
                                        0.0,
                                        leOMP->GetNumberOfFrameState(iChannel) != 0 ? leOMP->GetState(iChannel, t) : 0.0);
                                }
                            }else{
                                for(size_t t=0;t<leOMP->GetNumberOfFrameX(iChannel);t++){
                                    fprintf(fout, " %23.16f %23.16f %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n", 
                                        saveTimeStep * t, 
                                        leOMP->GetX(iChannel, t),
                                        leOMP->GetV(iChannel, t),
                                        leOMP->GetF(iChannel, t),
                                        leOMP->GetBarrier(iChannel, t),
                                        leOMP->GetPE(iChannel, t),
                                        leOMP->GetKE(iChannel, t),
                                        leOMP->GetState(iChannel, t));
                                }
                            }
                            if( bWaitingTimeDistribution ){
                                fprintf(fout, "# event tick\n");
                                for(size_t t=0;t<leOMP->GetBarrierTick(iChannel).size();t++){
                                    // internal ticks are in unit of timestep
                                    fprintf(fout, "# %10zd %23.16E\n", t, leOMP->GetBarrierTick(iChannel)[t] * timestep);
                                }
                            }
                            fclose(fout);
                        }

                        if( bFirstPassageTime ){
                            firstPassageTime[iMPI].push_back(leOMP->GetFirstPassageTimeAverage(iChannel));
                            firstPassageTime2[iMPI].push_back(leOMP->GetFirstPassageTimeAverage2(iChannel));
                            leOMP->GetInitialBarrierHeightHistogram(iChannel, &histInitialBarrierHeight[iMPI]);
                            leOMP->GetInstantBarrierHeightHistogram(iChannel, &histInstantBarrierHeight[iMPI]);
                        }

                        if( bTransitionPathTime ){
                            transitionPathTime[iMPI].push_back(leOMP->GetTransitionPathTimeAverage(iChannel));
                            transitionPathTime2[iMPI].push_back(leOMP->GetTransitionPathTimeAverage2(iChannel));
                            transitionPathTime2All.insert(transitionPathTime2All.end(), leOMP->GetTransitionPathTimeAll2(iChannel).begin(), leOMP->GetTransitionPathTimeAll2(iChannel).end());
                        }

                        if( bFirstPassageTimeCG ){
                            for(int jjj=0;jjj<firstPassageTimeCGSize.size();jjj++){
                                firstPassageTimeCG2[iMPI][jjj].push_back(leOMP->GetFirstPassageTimeCGAverage2(iChannel, jjj));
                            }
                        }

                        double stateCount0[2] = { 0, 0 };
                        stateCount0[0] = leOMP->GetStateCount(iChannel, 0);
                        stateCount0[1] = leOMP->GetStateCount(iChannel, 1);

                        stateCount[0][iMPI] += stateCount0[0];
                        stateCount[1][iMPI] += stateCount0[1];
                    }
                    timer2.StopToEntry("analysis");
                }

                loopTimer.Update(i, ntrial);
                fprintf(stderr, "   Done %10d/%10d %s\n", i + 1, ntrial, loopTimer.GetPChar());

                // synchronize
                // m_MPI.Barrier();
            }

            // m_MPI.Barrier();
            timer2.Start();

            if( bFirstPassageTime ){
                std::vector<double> firstPassageTimeCombined;
                std::vector<double> firstPassageTimeCombined2;
                for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                    firstPassageTimeCombined.insert(firstPassageTimeCombined.end(), firstPassageTime[iMPI].begin(), firstPassageTime[iMPI].end());
                    firstPassageTimeCombined2.insert(firstPassageTimeCombined2.end(), firstPassageTime2[iMPI].begin(), firstPassageTime2[iMPI].end());
                }
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<histInitialBarrierHeight[0].count.size();j++){
                        histInitialBarrierHeight[0].count[j] += histInitialBarrierHeight[iMPI].count[j];
                        histInstantBarrierHeight[0].count[j] += histInstantBarrierHeight[iMPI].count[j];
                    }
                }

                {
                    text filename = outputBase + "_firstPassageTime.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<firstPassageTimeCombined.size();index++){
                        writer.WriteDouble(firstPassageTimeCombined[index]);
                    }
                    writer.Close();
                }
                {
                    text filename = outputBase + "_firstPassageTime2.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<firstPassageTimeCombined2.size();index++){
                        writer.WriteDouble(firstPassageTimeCombined2[index]);
                    }
                    writer.Close();
                }

                histInitialBarrierHeight[0].NormalizePercent();
                histInstantBarrierHeight[0].NormalizePercent();

                {
                    {
                        text filename = outputBase + "_firstPassageTime.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&firstPassageTimeCombined));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&firstPassageTimeCombined));
                        fprintf(fout, "# count   %23zd\n", firstPassageTimeCombined.size());
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_firstPassageTime2.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&firstPassageTimeCombined2));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&firstPassageTimeCombined2));
                        fprintf(fout, "# count   %23zd\n", firstPassageTimeCombined2.size());
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_initialBarrierHeightHistogram.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        histInitialBarrierHeight[0].Print(fout);
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_instantBarrierHeightHistogram.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        histInstantBarrierHeight[0].Print(fout);
                        fclose(fout);
                    }
                }
            }

            {
                // state count
                double stateCountCombined[2] = { 0, 0 };
                for(int i=0;i<numberOfMPIThreads;i++){
                    stateCountCombined[0] += stateCount[0][i];
                    stateCountCombined[1] += stateCount[1][i];
                }

                {
                    text filename = outputBase + "_stateCount.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE* fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    fprintf(fout, "# state 0 %23.0f %13.7f\n", stateCountCombined[0], stateCountCombined[0] / (stateCountCombined[0] + stateCountCombined[1]));
                    fprintf(fout, "# state 1 %23.0f %13.7f\n", stateCountCombined[1], stateCountCombined[1] / (stateCountCombined[0] + stateCountCombined[1]));
                    fclose(fout);
                }
            }

            if( bTransitionPathTime ){
                std::vector<double> transitionPathTimeCombined;
                std::vector<double> transitionPathTimeCombined2;
                for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                    transitionPathTimeCombined.insert(transitionPathTimeCombined.begin(), transitionPathTime[iMPI].begin(), transitionPathTime[iMPI].end());
                    transitionPathTimeCombined2.insert(transitionPathTimeCombined2.begin(), transitionPathTime2[iMPI].begin(), transitionPathTime2[iMPI].end());
                }

                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<histInitialBarrierHeight[0].count.size();j++){
                        histInitialBarrierHeight[0].count[j] += histInitialBarrierHeight[iMPI].count[j];
                        histInstantBarrierHeight[0].count[j] += histInstantBarrierHeight[iMPI].count[j];
                    }
                }

                {
                    text filename = outputBase + "_transitionPathTime.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<transitionPathTimeCombined.size();index++){
                        writer.WriteDouble(transitionPathTimeCombined[index]);
                    }
                    writer.Close();
                }
                {
                    text filename = outputBase + "_transitionPathTime2.dat";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    StreamerWriter writer;
                    writer.Open(filename, std::ios::binary);
                    for(int index=0;index<transitionPathTimeCombined2.size();index++){
                        writer.WriteDouble(transitionPathTimeCombined2[index]);
                    }
                    writer.Close();
                }

                {
                    {
                        text filename = outputBase + "_transitionPathTime.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&transitionPathTimeCombined));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&transitionPathTimeCombined));
                        fprintf(fout, "# count   %23zd\n", transitionPathTimeCombined.size());
                        fclose(fout);
                    }
                    {
                        text filename = outputBase + "_transitionPathTime2.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&transitionPathTimeCombined2));
                        fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&transitionPathTimeCombined2));
                        fprintf(fout, "# count   %23zd\n", transitionPathTimeCombined2.size());
                        fclose(fout);
                    }
                }
                {
                    text filename = outputBase + "_transitionPathTimeHistogram.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE* fout = fopen(filename.c_str(), "w+");
                    double maximum = STATISTICS::MAX(&transitionPathTime2All);
                    HISTOGRAM::HISTOGRAM hist;
                    hist.Allocate(0.0, maximum * 1.05, 200);
                    hist.HistogramVector(&transitionPathTime2All);
                    hist.NormalizePercent();
                    PrintParameter(fout);
                    hist.Print(fout);
                    fclose(fout);
                }
            }

            if( bFirstPassageTimeCG ){
                std::vector<std::vector<double > > firstPassageTimeCGCombined2(firstPassageTimeCGSize.size());
                for(int j=0;j<firstPassageTimeCGSize.size();j++){
                    for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                        firstPassageTimeCGCombined2[j].insert(firstPassageTimeCGCombined2[j].begin(), firstPassageTimeCG2[iMPI][j].begin(), firstPassageTimeCG2[iMPI][j].end());
                    }
                }

                {
                    for(int j=0;j<firstPassageTimeCGSize.size();j++){
                        char buffer[256];
                        sprintf(buffer, "%s_firstPassageTimeCG2_%04d_%04d.dat", 
                            outputBase.c_str(), j, firstPassageTimeCGSize[j]);
                        text filename = buffer;
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        StreamerWriter writer;
                        writer.Open(filename, std::ios::binary);
                        for(int index=0;index<firstPassageTimeCGCombined2[j].size();index++){
                            writer.WriteDouble(firstPassageTimeCGCombined2[j][index]);
                        }
                        writer.Close();
                    }
                }

                {
                    {
                        text filename = outputBase + "_firstPassageTimeCG2.log";
                        MOLUTILITY::findNextBackUpFileName(filename.string());
                        FILE* fout = fopen(filename.c_str(), "w+");
                        PrintParameter(fout);
                        for(int j=0;j<firstPassageTimeCGSize.size();j++){
                            fprintf(fout, "# %4d %4d\n", j, firstPassageTimeCGSize[j]);
                            fprintf(fout, "# average %23.16E\n", STATISTICS::AVERAGE(&firstPassageTimeCGCombined2[j]));
                            fprintf(fout, "# sd      %23.16E\n", STATISTICS::RMSD(&firstPassageTimeCGCombined2[j]));
                            fprintf(fout, "# count   %23zd\n", firstPassageTimeCGCombined2[j].size());
                        }
                        fclose(fout);
                    }
                }
            }

            // v0, eta0
            if( bGetV0Distribution ){
                // mpi sum
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<histV0[0].count.size();j++){
                        histV0[0].count[j] += histV0[iMPI].count[j];
                    }
                }

                {
                    histV0[0].NormalizePercent();

                    text filename = outputBase + "_V0Distribution.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    histV0[0].Print(fout);
                    fclose(fout);
                }
            }

            if( bGetEta0Distribution ){
                // mpi sum
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<histEta0[0].count.size();j++){
                        histEta0[0].count[j] += histEta0[iMPI].count[j];
                    }
                }

                {
                    histEta0[0].NormalizePercent();

                    text filename = outputBase + "_Eta0Distribution.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    histEta0[0].Print(fout);
                    fclose(fout);
                }
            }

            if( bGetState0Distribution ){
                // mpi sum
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<histState0[0].count.size();j++){
                        histState0[0].count[j] += histState0[iMPI].count[j];
                    }
                }

                {
                    histState0[0].NormalizePercent();

                    text filename = outputBase + "_State0Distribution.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    histState0[0].Print(fout);
                    fclose(fout);
                }
            }

            // now this is for eta
            if( writeEtaResult ){
                for(int i=0;i<probTime.size();i++){
                    for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                        for(int j=0;j<hist[iMPI].m_Data.size();j++){
                            hist[0].m_Data[i][j] += hist[iMPI].m_Data[i][j];
                        }
                    }
                }

                {
                    hist[0].NormalizePercent();

                    text filename = outputBase + "_EtaPDF.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int i=0;i<probTime.size();i++){
                        fprintf(fout, "# prob t = %13.7f %10zd\n", probTime[i], probTimeStep[i]);
                    }
                    hist[0].Print(fout);
                    fclose(fout);
                }

                // average corrEta
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<corrEta[iMPI].size();j++){
                        corrEta[0][j] += corrEta[iMPI][j];
                    }
                }
                {
                    for(size_t t=0;t<interval.size();t++){
                        corrEta[0][t] /= double(nrun);
                    }

                    text filename = outputBase + "_EtaCorrelation.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int t=0;t<interval.size();t++){
                        fprintf(fout, "%13.7f %13.7f\n", saveTimeStep * interval[t], corrEta[0][t]);
                    }
                    fclose(fout);
                }
            }


            // x correlation
            if( bXCorr ) {
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<corr[iMPI].size();j++){
                        corr[0][j] += corr[iMPI][j];
                    }
                }

                {
                    for(size_t t=0;t<interval.size();t++){
                        corr[0][t] /= double(nrun);
                    }

                    text filename = outputBase + "_correlation.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int t=0;t<interval.size();t++){
                        fprintf(fout, "%13.7f %13.7f\n", saveTimeStep * interval[t], corr[0][t]);
                    }
                    fclose(fout);
                }
            }

            // msd
            {
                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<msd[iMPI].size();j++){
                        msd[0][j] += msd[iMPI][j];
                    }
                }

                {
                    for(size_t t=0;t<interval.size();t++){
                        msd[0][t] /= double(nrun);
                    }

                    text filename = outputBase + "_msd.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    for(int t=0;t<interval.size();t++){
                        fprintf(fout, "%14.7E %14.7E\n", saveTimeStep * interval[t], msd[0][t]);
                    }
                    fclose(fout);
                }
            }

            // crossing event
            if( bWaitingTimeDistribution ){
                std::vector<double> waitingTimeAvgS(numberOfMPIThreads);
                std::vector<double> waitingTimeAvg2S(numberOfMPIThreads);
                std::vector<double> waitingTimeMaxS(numberOfMPIThreads);
                std::vector<double> waitingTimeListCountS(numberOfMPIThreads);

                for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                    waitingTimeAvgS[iMPI] = STATISTICS::SUM(&waitingTimeList[iMPI]);
                    waitingTimeAvg2S[iMPI] = STATISTICS::SUM2(&waitingTimeList[iMPI]);
                    waitingTimeMaxS[iMPI] = STATISTICS::MAX(&waitingTimeList[iMPI]);
                    waitingTimeListCountS[iMPI] = waitingTimeList[iMPI].size();
                }

                double waitingTimeAvg = STATISTICS::SUM(&waitingTimeAvgS);
                double waitingTimeAvg2 = STATISTICS::SUM(&waitingTimeAvg2S);
                double waitingTimeListCount = STATISTICS::SUM(&waitingTimeListCountS);
                waitingTimeAvg /= waitingTimeListCount;
                waitingTimeAvg2 /= waitingTimeListCount;
            
                double waitingTimeMax = STATISTICS::MAX(&waitingTimeMaxS);

                double R = (waitingTimeAvg2 - waitingTimeAvg * waitingTimeAvg) / waitingTimeAvg / waitingTimeAvg;

                std::vector<HISTOGRAM::HISTOGRAM> histWaitingTime(numberOfMPIThreads);

                for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                    int nbin = 100;
                    if( waitingTimeListCount < 1000 ){
                        nbin = std::max<int>(10, (int)waitingTimeListCount / 100);
                    }
                    histWaitingTime[iMPI].Allocate(0.0, waitingTimeMax + 0.1, nbin);
                    histWaitingTime[iMPI].HistogramVector(&waitingTimeList[iMPI]);
                }

                for(int iMPI=1;iMPI<numberOfMPIThreads;iMPI++){
                    for(int j=0;j<histWaitingTime[iMPI].count.size();j++){
                        histWaitingTime[0].count[j] += histWaitingTime[iMPI].count[j];
                    }
                }

                {
                    histWaitingTime[0].NormalizePercent();

                    text filename = outputBase + "_waitingTime.log";
                    MOLUTILITY::findNextBackUpFileName(filename.string());
                    FILE *fout = fopen(filename.c_str(), "w+");
                    PrintParameter(fout);
                    fprintf(fout, "# R      = %13.7f\n", R);
                    fprintf(fout, "# Nwt    = %13zd\n", (size_t)waitingTimeListCount);
                    histWaitingTime[0].Print(fout);
                    fclose(fout);
                }
            }

            SAFE_DELETE(leOMP);
            timer2.StopToEntry("post");

            timer.StopToEntry("total");            

        };
    };

    // =================================================================================================
    // statistics
    // =================================================================================================
    class LangevinEquationRunnerMPIPostProcessStatistics{
    public:
        text outputBase = "langevinEquation";

        double x0;
        double v0;
        double gamma;
        double kb;
        double temperature;
        double mass[2];
        double timestep;
        size_t numberOfStep;

        size_t statisticsFreq;
        size_t saveFreq;

        bool bRandomizeVelocity;
        bool bRandomizeEta;
        bool bRandomizeState;
        bool bSameInitialBarrier;

        double periodicLength;
        double barrierHeight[2];
        int initialState;
        double slopeHalfLength;
        double xc;
        double yc;

        double q;
        double tau;
        double eta0;

        bool bWaitingTimeDistribution;
        bool bFirstPassageTime;
        bool bTransitionPathTime;
        bool bFirstPassageTimeCG;
        std::vector<int> firstPassageTimeCGSize;

        bool bReadRestart;
        // bool bWriteRestart;
        std::vector<text> sInputRestart[2];
        // text sOutputRestart;

        int ntrial;
        size_t numberOfStepPrevious;

        bool writeEtaResult;
        double probHMin;
        double probHMax;
        double probHBin;
        std::vector<double> probTime;

        bool writeX;
        int samplePerOrder;
        int randomNumberSeedForSeed;

        bool bXCorr;
        bool bTimer;

        bool bGetV0Distribution;
        bool bGetEta0Distribution;
        bool bGetState0Distribution;

        bool bOverdamped;

        HISTOGRAM::HistogramSetting m_BarrierHistogramSetting;
        //
        // MPIKERNEL::MPIKernel m_MPI;
        int numberOfMPIThreads;

        bool bEqualProbability;
        double equilibriumConstant;
    protected:
        int nrunPerNode;
        int nrun;
        int randomNumberSeedUsed;
        double numberOfTrajectorySet;
        double memoryUsage;
        double saveTimeStep;
    public:
        LangevinEquationRunnerMPIPostProcessStatistics(){
            outputBase = "langevinEquation";
            x0 = 0.0;
            v0 = 0.0;
            gamma = 1.0;
            kb = 1.0;
            temperature = 1.0;
            mass[0] = 1.0;
            mass[1] = 1.0;
            timestep = 0.1;
            numberOfStep = 10000;

            statisticsFreq = 1;
            saveFreq = 1;

            bRandomizeVelocity = false;
            bRandomizeEta = false;
            bRandomizeState = true;

            periodicLength = 5.0;
            barrierHeight[0] = 2.0;
            barrierHeight[1] = 6.0;
            initialState = 0;
            slopeHalfLength = 1.0;
            xc = 2.5;
            yc = 2.0;

            q = 0.2;
            tau = 0.1;
            eta0 = 0.0;
            
            bWaitingTimeDistribution = false;
            bFirstPassageTime = false;
            bTransitionPathTime = false;
            bFirstPassageTimeCG = false;

            bReadRestart = false;
            // bWriteRestart = false;
            // sInputRestart[0] = "inputRestart";
            // sInputRestart[1] = "inputRestart";
            // sOutputRestart = "outputRestart";

            ntrial = 100000;
            numberOfStepPrevious = 0;

            writeEtaResult = false;
            probHMin = -100.0;
            probHMax = 100.0;
            probHBin = 0.1;
            probTime = {0.1, 1.0, 10.0, 100.0};

            writeX = false;
            samplePerOrder = 20;
            randomNumberSeedForSeed = -1;

            bXCorr = false;
            bTimer = false;

            bGetV0Distribution = false;
            bGetEta0Distribution = false;
            bGetState0Distribution = false;

            bOverdamped = false;

            m_BarrierHistogramSetting.Set(0.0, 10.0, 100);

            // m_MPI.Initialize(m_World);

            numberOfMPIThreads = 1;

            nrunPerNode = 1;
            nrun = 1;
            randomNumberSeedUsed = -1;
            numberOfTrajectorySet = 0;
            memoryUsage = 0;
            saveTimeStep = 0.0;

            bEqualProbability = true;
            equilibriumConstant = 1.0;
        }
    public:
        void PrintParameter(FILE* fout = stderr){
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Overdamped           = %13s\n", text(bOverdamped).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Number of node       = %13d\n", numberOfMPIThreads);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Time step            = %13.7f\n", timestep);
            fprintf(fout, "# Save frequency       = %13zd\n", saveFreq);
            fprintf(fout, "# Save time step       = %13.7f\n", saveTimeStep);
            fprintf(fout, "# Statistics freq      = %13zd\n", statisticsFreq);
            fprintf(fout, "# Number of step (pre) = %13zd\n", numberOfStepPrevious);
            fprintf(fout, "# Number of step       = %13zd\n", numberOfStep);
            fprintf(fout, "# Simulation time      = %13.6E\n", timestep * numberOfStep);
            fprintf(fout, "# Number of trial      = %13d\n", ntrial);
            fprintf(fout, "# Number of run/trial  = %13d\n", SIMDWIDTH);
            fprintf(fout, "# Number of run        = %13d\n", nrun);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Barrier height [0]   = %13.7f\n", barrierHeight[0]);
            fprintf(fout, "# Barrier height [1]   = %13.7f\n", barrierHeight[1]);
            fprintf(fout, "# Xc                   = %13.7f\n", xc);
            fprintf(fout, "# Yc                   = %13.7f\n", yc);
            fprintf(fout, "# Barrier periodicity  = %13.7f\n", periodicLength);
            fprintf(fout, "# Slope half length    = %13.7f\n", slopeHalfLength);
            fprintf(fout, "# Equal probability    = %13s\n", text(bEqualProbability).c_str());
            fprintf(fout, "# Equilibrium constant = %13.7E\n", equilibriumConstant);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Mass [0]             = %13.7f\n", mass[0]);
            fprintf(fout, "# Mass [1]             = %13.7f\n", mass[1]);
            fprintf(fout, "# Gamma                = %13.7f\n", gamma);
            fprintf(fout, "# Temperature          = %13.7f\n", temperature);
            fprintf(fout, "# kB                   = %13.7f\n", kb);
            fprintf(fout, "# Eta variance (Q)     = %13.7f\n", q);
            fprintf(fout, "# Eta tau              = %13.7f\n", tau);
            fprintf(fout, "# Sample per order     = %13d\n", samplePerOrder);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# X0                   = %13.7f\n", x0);
            if( bRandomizeVelocity ){
                fprintf(fout, "# V0                   = %13s\n", "Random");
            }else{
                fprintf(fout, "# V0                   = %13.7f\n", v0);
            }
            if( bRandomizeEta ){
                fprintf(fout, "# Eta0                 = %13s\n", "Random");
            }else{
                fprintf(fout, "# Eta0                 = %13.7f\n", eta0);
            }
            if( bRandomizeState ){
                fprintf(fout, "# State0               = %13s\n", "Random");
            }else{
                fprintf(fout, "# State0               = %13d\n", initialState);
            }
            fprintf(fout, "# Same initial barrier = %13s\n", text(bSameInitialBarrier).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Random number seed0  = %13d\n", randomNumberSeedForSeed);
            fprintf(fout, "# Random number seed   = %13d\n", randomNumberSeedUsed);
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Number of traj(each) = %13.0f\n", numberOfTrajectorySet);
            fprintf(fout, "# Memory usage (each)  = %13s\n", text::FromByte2Readable(memoryUsage).c_str());
            fprintf(fout, "# Memory usage (total) = %13s\n", text::FromByte2Readable(memoryUsage * numberOfMPIThreads).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# MSD                  = %13s\n", "Mandatory");
            fprintf(fout, "# Waiting time         = %13s\n", text(bWaitingTimeDistribution).c_str());
            fprintf(fout, "# Write trajectories   = %13s\n", text(writeX).c_str());
            fprintf(fout, "# V0 Distribution      = %13s\n", text(bGetV0Distribution).c_str());
            fprintf(fout, "# Eta0 Distribution    = %13s\n", text(bGetEta0Distribution).c_str());
            fprintf(fout, "# State0 Distribution  = %13s\n", text(bGetState0Distribution).c_str());
            fprintf(fout, "# Eta result           = %13s\n", text(writeEtaResult).c_str());
            fprintf(fout, "# First passage time   = %13s\n", text(bFirstPassageTime).c_str());
            if( bFirstPassageTime ){
                fprintf(fout, "# Barrier hist hmin    = %13.7f\n", m_BarrierHistogramSetting.m_HistHMin);
                fprintf(fout, "# Barrier hist hmax    = %13.7f\n", m_BarrierHistogramSetting.m_HistHMax);
                fprintf(fout, "# Barrier hist nbin    = %13d\n", m_BarrierHistogramSetting.m_HistNBin);
            }
            fprintf(fout, "# Transtion path time  = %13s\n", text(bTransitionPathTime).c_str());
            fprintf(fout, "# First passage time CG= %13s\n", text(bFirstPassageTimeCG).c_str());
            if( bFirstPassageTimeCG ){
                for(int j=0;j<firstPassageTimeCGSize.size();j++){
                    fprintf(fout, "# FFT CG size %4d     = %13d\n", j, firstPassageTimeCGSize[j]);
                }
            }
            fprintf(fout, "# X correlation        = %13s\n", text(bXCorr).c_str());
            fprintf(fout, "# ==============================================\n");
            fprintf(fout, "# Read restart         = %13s\n", text(bReadRestart).c_str());
            fprintf(fout, "# Write restart        = %13s\n", "false");
            for(int i=0;i<sInputRestart[0].size();i++){
                fprintf(fout, "# Input restart (base0)= %13s\n", sInputRestart[0][i].c_str());
            }
            for(int i=0;i<sInputRestart[1].size();i++){
                fprintf(fout, "# Input restart (base1)= %13s\n", sInputRestart[1][i].c_str());
            }
            fprintf(fout, "# Output restart (base)= %13s\n", "none");
            fprintf(fout, "# ==============================================\n");
        }
        void Setup(
            text outputBase = "langevinEquation", 
            double x0 = 0.0, double v0 = 0.0, 
            double gamma = 1.0, double kb = 1.0, double temperature = 1.0, double mass0 = 1.0, double mass1 = 1.0, 
            double timestep = 0.1, size_t numberOfStep = 10000, 
            size_t statisticsFreq = 100, size_t saveFreq = 1, 
            bool bRandomizeVelocity = false, bool bRandomizeEta = false, bool bRandomizeState = true, 
            bool bSameInitialBarrier = false,
            double periodicLength = 5.0, 
            double barrierHeightBL = 1.0, double barrierHeightBH = 1.0, 
            int initialState = 0, 
            double slopeHalfLengthB = 1.0, 
            double q = 0.2, double tau = 0.1, double eta0 = 0.0,
            bool bWaitingTimeDistribution = false, double bFirstPassageTime = false,
            bool bReadRestart = false, // bool bWriteRestart = false, 
            std::vector<text> sInputRestart0 = {"inputRestart"}, // text sOutputRestart = "outputRestart", 
            std::vector<text> sInputRestart1 = {"inputRestart"}, // text sOutputRestart = "outputRestart", 
            int ntrial = 100000, size_t numberOfStepPrevious = 0, 
            bool writeEtaResult = false, double probHMin = -100.0, double probHMax = 100.0, 
            double probHBin = 0.1, std::vector<double> probTime = {0.1, 1.0, 10.0, 100.0},
            bool writeX = false, 
            int samplePerOrder = 20, int randomNumberSeedForSeed = -1,
            bool bXCorr = false, bool bTimer = false,
            bool bGetV0Distribution = false, bool bGetEta0Distribution = false, bool bGetState0Distribution = false,
            bool bOverdamped = false, 
            double barrierHistHMin = 0.0, double barrierHistHMax = 10.0, int barrierHistNBin = 100,
            bool bTransitionPathTime = false, bool bFirstPassageTimeCG = false, std::vector<int> firstPassageTimeCGSize = {},
            bool bEqualProbability = true, double equilibriumConstant = 1.0, int numberOfMPIThreads = 1){
            // code
            this->outputBase = outputBase;
            this->x0 = x0;
            this->v0 = v0;
            this->gamma = gamma;
            this->kb = kb;
            this->temperature = temperature;
            this->mass[0] = mass0;
            this->mass[1] = mass1;
            this->timestep = timestep;
            this->numberOfStep = numberOfStep;

            this->statisticsFreq = statisticsFreq;
            this->saveFreq = saveFreq;

            this->bRandomizeVelocity = bRandomizeVelocity;
            this->bRandomizeEta = bRandomizeEta;
            this->bRandomizeState = bRandomizeState;
            this->bSameInitialBarrier = bSameInitialBarrier;

            this->periodicLength = periodicLength;
            this->barrierHeight[0] = barrierHeightBL;
            this->barrierHeight[1] = barrierHeightBH;
            this->initialState = initialState;
            this->slopeHalfLength = slopeHalfLengthB;
            
            this->q = q;
            this->tau = tau;
            this->eta0 = eta0;

            this->bWaitingTimeDistribution = bWaitingTimeDistribution;
            this->bFirstPassageTime = bFirstPassageTime;
            this->bTransitionPathTime = bTransitionPathTime;
            this->bFirstPassageTimeCG = bFirstPassageTimeCG;
            this->firstPassageTimeCGSize = firstPassageTimeCGSize;

            this->bReadRestart = bReadRestart;
            // this->bWriteRestart = bWriteRestart;
            this->sInputRestart[0] = sInputRestart0;
            this->sInputRestart[1] = sInputRestart1;
            // this->sOutputRestart = sOutputRestart;

            this->ntrial = ntrial;
            this->numberOfStepPrevious = numberOfStepPrevious;

            this->writeEtaResult = writeEtaResult;
            this->probHMin = probHMin;
            this->probHMax = probHMax;
            this->probHBin = probHBin;
            this->probTime = probTime;

            this->writeX = writeX;
            this->samplePerOrder = samplePerOrder;
            this->randomNumberSeedForSeed = randomNumberSeedForSeed;

            this->bXCorr = bXCorr;
            this->bTimer = bTimer;

            this->bGetV0Distribution = bGetV0Distribution;
            this->bGetEta0Distribution = bGetEta0Distribution;
            this->bGetState0Distribution = bGetState0Distribution;

            this->bOverdamped = bOverdamped;

            this->m_BarrierHistogramSetting.Set(barrierHistHMin, barrierHistHMax, barrierHistNBin);

            this->bEqualProbability = bEqualProbability;
            this->equilibriumConstant = equilibriumConstant;

            this->numberOfMPIThreads = numberOfMPIThreads;
        }
    protected:
        std::vector<std::vector<int>> generateBootstrapVectors(std::mt19937& rng, int nset, int nsub, int nsample) {
            if (nsample <= 0 || nsub <= 0 || nset <= 0) {
                throw std::invalid_argument("All arguments must be positive");
            }

            std::uniform_int_distribution<int> dist(0, nsample - 1);
            std::unordered_set<std::string> seen;
            std::vector<std::vector<int>> result;

            while ((int)result.size() < nset) {
                std::vector<int> vec(nsub);
                for (int& v : vec) {
                    v = dist(rng);  // sample with replacement
                }

                // Optionally ensure each full vector is unique across result
                std::string key = vecToString(vec);
                if (seen.count(key) == 0) {
                    seen.insert(key);
                    result.push_back(std::move(vec));
                }
            }

            return result;
        }
        std::string vecToString(const std::vector<int>& vec) {
            std::ostringstream oss;
            for (int v : vec) oss << v << ',';
            return oss.str();  // Used as key for uniqueness
        }
    public:
        void Run(double fraction, long nset, int ncore, double eps){
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::mt19937 engine(seed);

            std::vector<std::vector<double> > msdValue[2];
            std::vector<double> msdTime[2];
            for(int k=0;k<2;k++){
                fprintf(stderr, "Read for mass %14.7f\n", mass[k]);
                for(int i=0;i<sInputRestart[k].size();i++){
                    auto [msdTimeRead, msdValueRead] = ParseMSDFromRestartFile(k, i);
                    if( i == 0 ){
                        msdTime[k] = msdTimeRead;
                    }
                    msdValue[k].insert(msdValue[k].end(), msdValueRead.begin(), msdValueRead.end());
                }
            }
            
            std::vector<double> msdValueAvg[2];
            msdValueAvg[0].resize(msdTime[0].size());
            msdValueAvg[1].resize(msdTime[1].size());

            bool bFitted[2] = { false, false };
            double slope[2] = { 0.0, 0.0 };
            double intercept[2] = { 0.0, 0.0 };
            double correl[2] = { 0.0, 0.0 };
            int leftBound[2] = { -1, -1 };
            int rightBound[2] = { -1, -1 };
            double epsTest[2] = { eps, eps };

            for(int i=0;i<2;i++){
                for(int j=0;j<msdValue[i].size();j++){
                    #pragma omp parallel for num_threads(ncore)
                    for(int k=0;k<msdValueAvg[i].size();k++){
                        msdValueAvg[i][k] += msdValue[i][j][k];
                    }
                }
                for(int k=0;k<msdValueAvg[i].size();k++){
                    msdValueAvg[i][k] /= double(msdValue[i].size());
                }

                int maxIteration = 1000;
                for(int j=0;j<maxIteration;j++){
                    std::tie(bFitted[i], slope[i], intercept[i], correl[i], leftBound[i], rightBound[i]) = 
                        MOLUTILITY::AutoLinearFitting(msdTime[i], msdValueAvg[i], true, epsTest[i], 3, "", 1);
                    if( bFitted[i] ){
                        break;
                    }else{
                        epsTest[i] = eps * (j + 2);
                    }
                }

                TextFileWriter writer;
                text filename = text::formatC("test_msd_%d.log", i);
                writer.Open(filename);
                for(int a=0;a<msdTime[i].size();a++){
                    writer.Write("%24.16E %24.16E\n", msdTime[i][a], msdValueAvg[i][a]);
                }
                writer.Close();
            }

            long nsample = msdValue[0].size();
            long nsub = std::round(double(nsample) * fraction);
            std::vector<std::vector<int> > randomList = generateBootstrapVectors(engine, nset, nsub, nsample);

            // for(int i=0;i<nset;i++){
            //     fprintf(stderr, "%3d : %s\n", i, vecToString(randomList[i]).c_str());
            // }

            fprintf(stderr, "NSample %10ld NSub %10ld NTrial %10ld\n", nsample, nsub, nset);
            fprintf(stderr, "Avg D   %24.16E %24.16E\n", slope[0] * 0.5, slope[1] * 0.5);
            fprintf(stderr, "Avg D/D %24.16E\n", (slope[0] * mass[0]) / (slope[1] * mass[1]));
            fprintf(stderr, "Eps     %24.16E %24.16E\n", epsTest[0], epsTest[1]);
            fprintf(stderr, "Left B  %24.16E %24.16E\n", msdTime[0][leftBound[0]], msdTime[1][leftBound[1]]);
            fprintf(stderr, "Right B %24.16E %24.16E\n", msdTime[0][rightBound[0]], msdTime[1][rightBound[1]]);
            fprintf(stderr, "Fitted  %24s %24s\n", text(bFitted[0]).c_str(), text(bFitted[1]).c_str());

            std::vector<int> indexer[2];
            indexer[0] = VectorAuxiliary::GenerateRange<int>(nsample);
            indexer[1] = VectorAuxiliary::GenerateRange<int>(nsample);

            std::shuffle(indexer[0].begin(), indexer[0].end(), engine);
            std::shuffle(indexer[1].begin(), indexer[1].end(), engine);

            std::vector<double> diffusioCoefficient[2];
            std::vector<double> diffusionRatio;
            diffusioCoefficient[0].resize(nset);
            diffusioCoefficient[1].resize(nset);
            diffusionRatio.resize(nset);

            std::vector<double> msdTimeFit[2];
            for(int k=0;k<2;k++){
                auto first = msdTime[k].begin();
                auto last = msdTime[k].begin();
                if( bFitted[k] ){
                    first += leftBound[k];
                    last += rightBound[k];
                }

                msdTimeFit[k].assign(first, last);
            }

            #pragma omp parallel for num_threads(ncore)
            for(int i=0;i<nset;i++){
                std::vector<double> msdValue0[2];
                msdValue0[0].resize(msdTime[0].size());
                msdValue0[1].resize(msdTime[1].size());

                for(int j=0;j<nsub;j++){
                    int index = randomList[i][j];
                    int index0 = indexer[0][index];
                    int index1 = indexer[1][index];
                    for(int k=0;k<msdValue0[0].size();k++){
                        msdValue0[0][k] += msdValue[0][index0][k];
                        msdValue0[1][k] += msdValue[1][index1][k];
                    }
                }

                for(int k=0;k<msdValue0[0].size();k++){
                    msdValue0[0][k] /= double(nsub);
                    msdValue0[1][k] /= double(nsub);
                }

                double slopeI[2];
                for(int k=0;k<2;k++){
                    auto first = msdValue0[k].begin();
                    auto last = msdValue0[k].begin();
                    if( bFitted[k] ){
                        first += leftBound[k];
                        last += rightBound[k];
                    }

                    std::vector<double> data(first, last);

                    double resultSlope, resultIntercept, resultCorrel;
                    MOLUTILITY::linearRegression(&msdTimeFit[k], &data, resultSlope, resultIntercept, resultCorrel);

                    slopeI[k] = resultSlope;
                }

                diffusioCoefficient[0][i] = slopeI[0] * 0.5;
                diffusioCoefficient[1][i] = slopeI[1] * 0.5;

                diffusionRatio[i] = (diffusioCoefficient[0][i] * mass[0]) / (diffusioCoefficient[1][i] * mass[1]);
            }

            double diffusionCoefficientSD[2];
            double diffusionCoefficientRatioSD;
            diffusionCoefficientSD[0] = STATISTICS::RMSD(&diffusioCoefficient[0]);
            diffusionCoefficientSD[1] = STATISTICS::RMSD(&diffusioCoefficient[1]);
            diffusionCoefficientRatioSD = STATISTICS::RMSD(&diffusionRatio);

            fprintf(stderr, "SD  D   %24.16E %24.16E\n", diffusionCoefficientSD[0], diffusionCoefficientSD[1]);
            fprintf(stderr, "SD D/D  %24.16E\n", diffusionCoefficientRatioSD);

        }
        std::tuple<std::vector<double>, std::vector<std::vector<double> > > ParseMSDFromRestartFile(int massID, int restartID){
            // code
            fprintf(stderr, "LangevinEquationMultipleRunner\n");

            size_t numberOfFrame = LangevinEquationSIMDMultiple::CalculateNumberOfFrame(numberOfStep + numberOfStepPrevious, saveFreq);
            if( numberOfFrame < 1000 ){
                fprintf(stderr, "Error: numberOfFrame < 1000\n");
                return std::make_tuple(std::vector<double>(), std::vector<std::vector<double> >());
            }

            nrunPerNode = ntrial * SIMDWIDTH;
            nrun = nrunPerNode * numberOfMPIThreads;

            // std::vector<int> randomNumberSeedList = m_MPI.RandomSeed(randomNumberSeedForSeed);

            LangevinEquationSIMDMultiple* leOMP = 0;
            {
                if( bOverdamped ){
                    try{
                        leOMP = new OverdampedLangevinEquationSIMDMultiple();
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", 0, e.what());
                        exit(0);
                    }
                }else{
                    try{
                        leOMP = new LangevinEquationSIMDMultiple();
                    }catch( std::bad_alloc& e ){
                        fprintf(stderr, "Error: fail to allocate langevin equation on node %10d %s\n", 0, e.what());
                        exit(0);
                    }
                }

                leOMP->m_X0 = x0;
                leOMP->m_V0 = v0;
                leOMP->m_Gamma0 = gamma;
                leOMP->m_kB0 = kb;
                leOMP->m_T0 = temperature;
                leOMP->m_M0 = mass[massID];
                leOMP->m_TimeStep0 = timestep;
                leOMP->m_NumberOfStep = numberOfStep;

                leOMP->m_StatisticsFreq = statisticsFreq;
                leOMP->m_SaveFreq = saveFreq;

                leOMP->m_SaveX = true;
                if( writeX ){
                    leOMP->m_SaveV = true;
                    leOMP->m_SaveF = true;
                    leOMP->m_SaveE = true;
                    leOMP->m_SaveState = true;
                    leOMP->m_SaveBarrier = true;
                }
                if( writeEtaResult ){
                    leOMP->m_SaveEta = true;
                }

                leOMP->m_RandomizeVelocity = bRandomizeVelocity;
                leOMP->m_RandomizeEta = bRandomizeEta;
                leOMP->m_RandomizeState0 = bRandomizeState;
                leOMP->m_bSameInitialBarrier = bSameInitialBarrier;

                leOMP->m_PeriodicLength0 = periodicLength;
                leOMP->m_BarrierHeight0[0] = barrierHeight[0];
                leOMP->m_BarrierHeight0[1] = barrierHeight[1];
                leOMP->m_SlopeHalfLength0 = slopeHalfLength;

                leOMP->m_InitialState0 = initialState;
                leOMP->m_Q0 = q;
                leOMP->m_Tau0 = tau;
                leOMP->m_Eta0 = eta0;

                leOMP->m_bWaitingTimeDistribution = bWaitingTimeDistribution;
                leOMP->m_bFirstPassageTime = bFirstPassageTime;
                leOMP->m_bTransitionPathTime = bTransitionPathTime;
                leOMP->m_bFirstPassageTimeCG = bFirstPassageTimeCG;
                leOMP->m_FirstPassageTimeCGSize = firstPassageTimeCGSize;
                
                leOMP->m_ReadRestartFile = bReadRestart;
                leOMP->m_WriteRestartFile = false;
                leOMP->m_InputRestartFileName = sInputRestart[massID][restartID];
                leOMP->m_OutputRestartFileName = "none";

                leOMP->m_bEqualProbability = bEqualProbability;
                leOMP->m_EquilibriumConstant = equilibriumConstant;
                
                leOMP->m_Rank = 0;

                // this depends on tau 
                // randomNumberSeedUsed = leOMP->Initialize(randomNumberSeed);
                randomNumberSeedUsed = -1;
                leOMP->Initialize(randomNumberSeedUsed);

                // if( m_MPI.IsMasterNode() && ntrial != 1 && bReadRestart ){
                //     fprintf(stderr, "Warning: Exact restart can only be produced of the first trial\n");
                //     fprintf(stderr, "         Because the random number from the same seed is consumed differently\n");
                // }
            }

            numberOfTrajectorySet = leOMP->CalculateNumberOfTrajectorySet();
            memoryUsage = leOMP->CalculateMemoryUsage();

            saveTimeStep = timestep * saveFreq;

            // if( m_MPI.IsMasterNode() ){
                PrintParameter();
            // }

            // int maxlag = maxLagTime / timestep;
            std::vector<size_t> interval = 
                MOLUTILITY::GenerateLog10Scale(
                    saveTimeStep, numberOfFrame, (size_t)1000, (size_t)samplePerOrder, true);
            if( interval.size() < 20 ){
                fprintf(stderr, "Error: interval is less than 20 (%zd)\n", interval.size());
                SAFE_DELETE(leOMP);

                MPI_Finalize();
                exit(0);
            }else if( interval.size() < 100 ){
                fprintf(stderr, "Warning: interval is less than 100 (%zd)\n", interval.size());
            }else{
                fprintf(stderr, "Interval = %10zd\n", interval.size());
            }

            std::vector<double> intervalDouble(interval.size());
            for(int i=0;i<interval.size();i++){
                intervalDouble[i] = saveTimeStep * interval[i];
            }

            // msd (mandatory)
            std::vector<std::vector<double> > msd;
            VectorAuxiliary::AllocateVector2D(&msd, ntrial * numberOfMPIThreads * SIMDWIDTH, interval.size());

            clsSystemTimeCPP timer;
            timer.AddEntry("total");

            clsSystemTimeCPP timer2;
            timer2.AddEntry("integrate");
            timer2.AddEntry("analysis");
            timer2.AddEntry("post");

            timer.Start();

            LoopTimer loopTimer;
            loopTimer.Start();

            for(int i=0;i<ntrial;i++){
                fprintf(stderr, "##################################################################\n");
                fprintf(stderr, "Trial %10d/%10d\n", i, ntrial);
                fprintf(stderr, "##################################################################\n");
                
                for(int iMPI=0;iMPI<numberOfMPIThreads;iMPI++){
                    leOMP->m_Rank = iMPI;
                    fprintf(stderr, "   Read          %10d\n", iMPI);
                    timer2.Start();
                    if( !leOMP->ReadRestartFileAll(i, false) ){
                        fprintf(stderr, "Error: can't open restart file\n");
                        SAFE_DELETE(leOMP);
                        exit(0);
                    }
                    timer2.StopToEntry("integrate");

                    fprintf(stderr, "   Post Analysis %10d\n", iMPI);
                    timer2.Start();

                    for(int iChannel=0;iChannel<SIMDWIDTH;iChannel++){
                        int iMSD = i * numberOfMPIThreads * SIMDWIDTH + iMPI * SIMDWIDTH + iChannel;
                        // msd
                        MeanSquareDisplacementBruteForce(
                            &leOMP->GetX(iChannel), &msd[iMSD], &interval);

                        for(int t=0;t<interval.size();t++){
                            if( isnan(msd[iMSD][t]) || isinf(msd[iMSD][t]) ){
                                fprintf(stderr, "Error: on node %10d, ill msd0\n", iMPI);
                            }
                        }
                    }
                    timer2.StopToEntry("analysis");
                }

                loopTimer.Update(i, ntrial);
                fprintf(stderr, "   Done %10d/%10d %s\n", i + 1, ntrial, loopTimer.GetPChar());

                // synchronize
                // m_MPI.Barrier();
            }

            // m_MPI.Barrier();
            // timer2.Start();
            SAFE_DELETE(leOMP);          

            return std::make_tuple(intervalDouble, msd);
        };
    };

};

#endif
