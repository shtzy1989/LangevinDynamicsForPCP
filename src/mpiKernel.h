#ifndef MD_MPI_MPI_KERNEL_H
#define MD_MPI_MPI_KERNEL_H

#include "array.h"
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
#include <mpi.h>

namespace MPIKERNEL{
    class MPIKernel{
    protected:
        MPI_Comm m_World;
        int m_NumberOfNode;
        int m_NodeIndex;
        std::vector<text> m_Argument;
    public:
        void Initialize(int* argc, char** argv[]);
        void Initialize(MPI_Comm world);
        size_t GetArgumentSize();
        text GetArgument(int index);
        void PrintInformation(FILE *fout = stderr);
        void Barrier();
        int GetNumberOfNode();
        int GetNodeIndex();
        bool IsMasterNode();
        int GetMasterNodeIndex();
        MPI_Comm GetWorld();
        bool BCast(void *buffer, int count, MPI_Datatype datatype);
        bool BCast(void *buffer, int count, MPI_Datatype datatype, int root);
        bool AllReduce(const void* buffer, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op);
        bool Reduce(const void* buffer, void* recvbuf, int count, MPI_Datatype datatype, int root, MPI_Op op);
        void CleanUp();
        bool IsNodeValid();
    public:
        void SumVectorDouble(std::vector<double> *data);
        void SumVectorInt(std::vector<int> *data);
        std::vector<int> RandomSeed(int value, int* instanceValue = 0);
        std::vector<double> CombineVectorDoubleToMaster(std::vector<double> *data);
        std::vector<int> CombineVectorIntToMaster(std::vector<int> *data);
        void SumHistogramLeftAignOMP(std::vector<double> *output, std::vector<std::vector<double> > *input);
        bool SumHistogramCount(HISTOGRAM::HISTOGRAM *hist, int ncore = 1);
        bool GetHistogramSetting(HISTOGRAM::HISTOGRAM *hist, HISTOGRAM::HistogramSetting *setting);
        int PrintFromNode(FILE* stream, int nodeIndex, const char* format, ...);
        int Print(FILE* stream, const char* format, ...);
        int PrintFromLocalNode(FILE* stream, const char* format, ...);
        void SendRecv(void* buffer, int count, int src, int dest, MPI_Datatype dataType);
        void Send(const void* sendBuffer, int count, int dest, MPI_Datatype dataType);
        void Recv(void* recvBuffer, int count, int src, MPI_Datatype dataType);
        void AllGather(void* sendBuffer, int count, void* recvBuffer, MPI_Datatype dataType);
        void Allgatherv(void* sendBuffer, int sendCount, MPI_Datatype sendType, void* recvBuffer, int recvCount[], int recvDisp[], MPI_Datatype recvType);
        void Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root);
        void Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root);
    public:
        static std::vector<int> GetDisplacement(std::vector<int>& count);
    };

    class MPISub : public MPIKernel{
    protected:
        MPIKernel* m_Parent;
    public:
        bool InitializeByNumber(MPIKernel* mpi, int number);
        bool Split(MPIKernel* mpi, int color);
        MPIKernel* GetParent();
        void CleanSub();
    };
}

#endif
