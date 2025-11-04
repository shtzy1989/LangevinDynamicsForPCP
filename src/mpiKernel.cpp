#include "mpiKernel.h"

namespace MPIKERNEL{
    void MPIKernel::Initialize(int* argc, char** argv[]){
        MPI_Init(argc, argv);
        for(int i=1;i<*argc;i++){
            m_Argument.push_back((*argv)[i]);
        }
        Initialize(MPI_COMM_WORLD);
    };
    void MPIKernel::Initialize(MPI_Comm world){
        m_World = world;
        if( IsNodeValid() ){
            MPI_Comm_size(m_World, &m_NumberOfNode);
            MPI_Comm_rank(m_World, &m_NodeIndex);
        }else{
            m_NumberOfNode = -1;
            m_NodeIndex = -1;
        }
    };
    size_t MPIKernel::GetArgumentSize(){
        return m_Argument.size();
    };
    text MPIKernel::GetArgument(int index){
        return m_Argument[index];
    };
    void MPIKernel::PrintInformation(FILE *fout){
        fprintf(fout, "Node %10d of %10d nodes\n", m_NodeIndex, m_NumberOfNode);
    };
    void MPIKernel::Barrier(){
        if( IsNodeValid() ){
            MPI_Barrier(GetWorld());
        }
    };
    int MPIKernel::GetNumberOfNode(){
        return m_NumberOfNode;
    };
    int MPIKernel::GetNodeIndex(){
        return m_NodeIndex;
    };
    bool MPIKernel::IsMasterNode(){
        return m_NodeIndex == 0;
    };
    int MPIKernel::GetMasterNodeIndex(){
        return 0;
    };
    MPI_Comm MPIKernel::GetWorld(){
        return m_World;
    };
    bool MPIKernel::BCast(void *buffer, int count, MPI_Datatype datatype){
        if( IsNodeValid() ){
            MPI_Bcast(buffer, count, datatype, GetMasterNodeIndex(), GetWorld());
        }
        return true;
    };
    bool MPIKernel::BCast(void *buffer, int count, MPI_Datatype datatype, int root){
        if( IsNodeValid() ){
            MPI_Bcast(buffer, count, datatype, root, GetWorld());
        }
        return true;
    };
    bool MPIKernel::AllReduce(const void* buffer, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op){
        if( IsNodeValid() ){
            MPI_Allreduce(buffer, recvbuf, count, datatype, op, GetWorld());
        }
        return true;
    }
    bool MPIKernel::Reduce(const void* buffer, void* recvbuf, int count, MPI_Datatype datatype, int root, MPI_Op op){
        if( IsNodeValid() ){
            MPI_Reduce(buffer, recvbuf, count, datatype, op, root, GetWorld());
        }
        return true;
    }
    void MPIKernel::CleanUp(){
        MPI_Finalize();
    };
    bool MPIKernel::IsNodeValid(){
        return m_World != MPI_COMM_NULL;
    }
    //
    void MPIKernel::SumVectorDouble(std::vector<double> *data){
        if( IsNodeValid() ){
            int size = data->size();
            std::vector<int> sizeBuffer(GetNumberOfNode());
            MPI_Allgather(&size, 1, MPI_INT, &sizeBuffer[0], 1, MPI_INT, GetWorld());
            for(int i=0;i<sizeBuffer.size();i++){
                if( sizeBuffer[i] != size ){
                    fprintf(stderr, "Error: vector size not consistent (%d) vs (%d)\n", size, sizeBuffer[i]);
                }
            }

            std::vector<double> buffer(data->size());
            for(int r=1;r<GetNumberOfNode();r++){
                if( IsMasterNode() ){
                    MPI_Status status;
                    MPI_Recv(&buffer[0], data->size(), MPI_DOUBLE, r, 0, GetWorld(), &status);
                    for(int i=0;i<data->size();i++){
                        (*data)[i] += buffer[i];
                    }
                }else if( r == GetNodeIndex() ){
                    memcpy(&buffer[0], &(*data)[0], sizeof(double) * data->size());
                    MPI_Send(&buffer[0], data->size(), MPI_DOUBLE, GetMasterNodeIndex(), 0, GetWorld());
                }
            }
        }
    }    
    void MPIKernel::SumVectorInt(std::vector<int> *data){
        if( IsNodeValid() ){
            std::vector<int> buffer(data->size());
            for(int r=1;r<GetNumberOfNode();r++){
                if( IsMasterNode() ){
                    MPI_Status status;
                    MPI_Recv(&buffer[0], data->size(), MPI_INT, r, 0, GetWorld(), &status);
                    for(int i=0;i<data->size();i++){
                        (*data)[i] += buffer[i];
                    }
                }else if( r == GetNodeIndex() ){
                    memcpy(&buffer[0], &(*data)[0], sizeof(int) * data->size());
                    MPI_Send(&buffer[0], data->size(), MPI_INT, GetMasterNodeIndex(), 0, GetWorld());
                }
            }
        }
    }
    std::vector<int> MPIKernel::RandomSeed(int value, int* instanceValue){
        std::vector<int> randomNumberSeedList;

        if( IsNodeValid() ){
            if( GetNodeIndex() == 0 ){
                int value0 = RANDOMNUMBER::SEED(value);
                if( instanceValue != 0 ) *instanceValue = value0;
                randomNumberSeedList = RANDOMNUMBER::RANDOM_NON_REPEAT_INT(GetNumberOfNode());
            }else{
                randomNumberSeedList.resize(GetNumberOfNode());
            }
            MPI_Bcast(&randomNumberSeedList[0], GetNumberOfNode(), MPI_INT, 0, GetWorld());
        }
        return randomNumberSeedList;
    }
    std::vector<double> MPIKernel::CombineVectorDoubleToMaster(std::vector<double> *data){
        std::vector<double> result;
        if( IsNodeValid() ){
            for(int i=0;i<GetNumberOfNode();i++){
                // send/recv size
                size_t buffer = data->size();
                if( i == GetNodeIndex() && !IsMasterNode() ){
                    // send
                    MPI_Send(&buffer, 1, MPI_UNSIGNED_LONG, GetMasterNodeIndex(), GetNodeIndex(), GetWorld());
                }else if( i != GetNodeIndex() && IsMasterNode() ){
                    // recv
                    MPI_Status status;
                    MPI_Recv(&buffer, 1, MPI_UNSIGNED_LONG, i, i, GetWorld(), &status);
                }
                // send/recv data
                if( i == GetNodeIndex() && !IsMasterNode() ){
                    // send
                    MPI_Send(&(*data)[0], buffer, MPI_DOUBLE, GetMasterNodeIndex(), GetNodeIndex(), GetWorld());
                }else if( i != GetNodeIndex() && IsMasterNode() ){
                    // recv
                    std::vector<double> buffer2;
                    MPI_Status status;
                    buffer2.resize(buffer);
                    MPI_Recv(&buffer2[0], buffer, MPI_DOUBLE, i, i, GetWorld(), &status);
                    for(int j=0;j<buffer2.size();j++){
                        result.push_back(buffer2[j]);
                    }
                }else if( IsMasterNode() && i == GetNodeIndex() ){
                    for(int j=0;j<data->size();j++){
                        result.push_back((*data)[j]);
                    }
                }
            }
        }
        return result;
    }
    std::vector<int> MPIKernel::CombineVectorIntToMaster(std::vector<int> *data){
        std::vector<int> result;
        if( IsNodeValid() ){
            for(int i=0;i<GetNumberOfNode();i++){
                // send/recv size
                size_t buffer = data->size();
                if( i == GetNodeIndex() && !IsMasterNode() ){
                    // send
                    MPI_Send(&buffer, 1, MPI_UNSIGNED_LONG, GetMasterNodeIndex(), GetNodeIndex(), GetWorld());
                }else if( i != GetNodeIndex() && IsMasterNode() ){
                    // recv
                    MPI_Status status;
                    MPI_Recv(&buffer, 1, MPI_UNSIGNED_LONG, i, i, GetWorld(), &status);
                }
                // send/recv data
                if( i == GetNodeIndex() && !IsMasterNode() ){
                    // send
                    MPI_Send(&(*data)[0], buffer, MPI_INT, GetMasterNodeIndex(), GetNodeIndex(), GetWorld());
                }else if( i != GetNodeIndex() && IsMasterNode() ){
                    // recv
                    std::vector<int> buffer2;
                    MPI_Status status;
                    buffer2.resize(buffer);
                    MPI_Recv(&buffer2[0], buffer, MPI_INT, i, i, GetWorld(), &status);
                    for(int j=0;j<buffer2.size();j++){
                        result.push_back(buffer2[j]);
                    }
                }else if( IsMasterNode() && i == GetNodeIndex() ){
                    for(int j=0;j<data->size();j++){
                        result.push_back((*data)[j]);
                    }
                }
            }
        }
        return result;
    }

    void MPIKernel::SumHistogramLeftAignOMP(std::vector<double> *output, std::vector<std::vector<double> > *input){
        if( IsNodeValid() ){
            // reduce OMP
            int ncore = input->size();
            for(int c=0;c<ncore;c++){
                if( output->size() < (*input)[c].size() ){
                    output->resize((*input)[c].size());
                }
                for(int h=0;h<(*input)[c].size();h++){
                    (*output)[h] += (*input)[c][h];
                }
            }
            // reduce MPI
            int size = output->size();
            int maxsize;
            AllReduce(&size, &maxsize, 1, MPI_INT, MPI_MAX);
            output->resize(maxsize);
            SumVectorDouble(output);
        }
    }

    bool MPIKernel::SumHistogramCount(HISTOGRAM::HISTOGRAM *hist, int ncore){
        if( IsNodeValid() ){
            double min0, max0;
            std::vector<double> step0(GetNumberOfNode());

            AllReduce(&hist->min, &min0, 1, MPI_DOUBLE, MPI_MIN);
            AllReduce(&hist->max, &max0, 1, MPI_DOUBLE, MPI_MAX);
            MPI_Allgather(&hist->step, 1, MPI_DOUBLE, &step0[0], 1, MPI_DOUBLE, GetWorld());

            int bOK = 0;
            for(int i=0;i<GetNumberOfNode();i++){
                if( step0[i] != hist->step ){
                    bOK++;
                }
            }

            int bAllOK = 0;
            AllReduce(&bOK, &bAllOK, 1, MPI_INT, MPI_SUM);

            if( bAllOK == 0 ){
                int size = (max0 - min0) / hist->step;
                std::vector<double> count(size);
                int index0 = (hist->min - min0) / hist->step;
                #pragma omp parallel for num_threads(ncore)
                for(int i=0;i<hist->size;i++){
                    count[i + index0] = hist->count[i];
                }
                SumVectorDouble(&count);

                if( IsMasterNode() ){
                    hist->min = min0;
                    hist->max = max0;
                    hist->size = size;
                    hist->count = count;
                    hist->percent.resize(count.size());
                }
                return true;
            }else{
                return false;
            }
        }else{
            return true;
        }
    }
    bool MPIKernel::GetHistogramSetting(HISTOGRAM::HISTOGRAM *hist, HISTOGRAM::HistogramSetting *setting){
        if( IsNodeValid() ){
            double min0, max0;
            std::vector<double> step0(GetNumberOfNode());

            AllReduce(&hist->min, &min0, 1, MPI_DOUBLE, MPI_MIN);
            AllReduce(&hist->max, &max0, 1, MPI_DOUBLE, MPI_MAX);
            MPI_Allgather(&hist->step, 1, MPI_DOUBLE, &step0[0], 1, MPI_DOUBLE, GetWorld());

            int bOK = true;
            for(int i=0;i<GetNumberOfNode();i++){
                if( step0[i] != hist->step ){
                    bOK++;
                }
            }

            int bAllOK = true;
            AllReduce(&bOK, &bAllOK, 1, MPI_INT, MPI_SUM);

            if( bAllOK ){
                setting->Set(min0, max0, hist->step);
                return true;
            }else{
                return false;
            }
        }else{
            return true;
        }
    }
    int MPIKernel::PrintFromNode(FILE* stream, int nodeIndex, const char* format, ...){
        if( IsNodeValid() ){
            if( GetNodeIndex() == nodeIndex ){
                va_list args;
                va_start(args, format);

                // You can add custom behavior here. For example, logging.

                int result = vfprintf(stream, format, args);

                va_end(args);

                // Additional behavior can be added here as well.

                return result;
            }else{
                return 0;
            }
        }else{
            return 0;
        }
    }
    int MPIKernel::Print(FILE* stream, const char* format, ...){
        if( IsNodeValid() ){
            if( IsMasterNode() ){
                va_list args;
                va_start(args, format);

                // You can add custom behavior here. For example, logging.

                int result = vfprintf(stream, format, args);

                va_end(args);

                // Additional behavior can be added here as well.

                return result;
            }else{
                return 0;
            }
        }else{
            return 0;
        }
    }
    int MPIKernel::PrintFromLocalNode(FILE* stream, const char* format, ...){
        if( IsNodeValid() ){
            va_list args;
            va_start(args, format);

            // You can add custom behavior here. For example, logging.

            int result = vfprintf(stream, format, args);

            va_end(args);

            // Additional behavior can be added here as well.

            return result;
        }else{
            return 0;
        }
    }
    void MPIKernel::SendRecv(void* buffer, int count, int src, int dest, MPI_Datatype dataType){
        if( IsNodeValid() ){
            if( GetNodeIndex() == src ){
                MPI_Send(buffer, count, dataType, dest, src, GetWorld());
            }else if( GetNodeIndex() == dest ){
                MPI_Status status;
                MPI_Recv(buffer, count, dataType, src, src, GetWorld(), &status);
            }
        }
    }
    void MPIKernel::Send(const void* sendBuffer, int count, int dest, MPI_Datatype dataType){
        if( IsNodeValid() ){
            MPI_Send(sendBuffer, count, dataType, dest, GetNodeIndex(), GetWorld());
        }
    }
    void MPIKernel::Recv(void* recvBuffer, int count, int src, MPI_Datatype dataType){
        if( IsNodeValid() ){
            MPI_Status status;
            MPI_Recv(recvBuffer, count, dataType, src, src, GetWorld(), &status);
        }
    }
    void MPIKernel::AllGather(void* sendBuffer, int count, void* recvBuffer, MPI_Datatype dataType){
        if( IsNodeValid() ){
            MPI_Allgather(sendBuffer, count, dataType, recvBuffer, count, dataType, GetWorld());
        }
    }
    void MPIKernel::Allgatherv(void* sendBuffer, int sendCount, MPI_Datatype sendType, void* recvBuffer, int recvCount[], int recvDisp[], MPI_Datatype recvType){
        MPI_Allgatherv(sendBuffer, sendCount, sendType, recvBuffer, recvCount, recvDisp, recvType, GetWorld());
    }
    void MPIKernel::Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
        void *recvbuf, int recvcount, MPI_Datatype recvtype, int root){
        if( IsNodeValid() ){
            MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, GetWorld());
        }
    }
    void MPIKernel::Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
        void *recvbuf, int recvcount, MPI_Datatype recvtype, int root){
        if( IsNodeValid() ){
            MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, GetWorld());
        }
    }
    std::vector<int> MPIKernel::GetDisplacement(std::vector<int>& count){
        std::vector<int> displs(count.size(), 0);
        for(int i=1;i<count.size();i++){
            displs[i] = displs[i - 1] + count[i - 1];
        }
        return displs;
    }

    bool MPISub::InitializeByNumber(MPIKernel* mpi, int number){
        if( number <= mpi->GetNumberOfNode() ){
            m_Parent = mpi;

            MPI_Group worldGroup;
            MPI_Comm_group(m_Parent->GetWorld(), &worldGroup);

            std::vector<int> ranks;
            for(int i=0;i<number;i++){
                ranks.push_back(i);
            }

            MPI_Group subGroup;
            MPI_Group_incl(worldGroup, ranks.size(), &ranks[0], &subGroup);

            MPI_Comm subComm;
            MPI_Comm_create_group(m_Parent->GetWorld(), subGroup, 0, &subComm);

            Initialize(subComm);

            MPI_Group_free(&worldGroup);
            MPI_Group_free(&subGroup);

            return true;
        }else{
            return false;
        }
    }
    bool MPISub::Split(MPIKernel* mpi, int color){
        m_Parent = mpi;

        MPI_Comm subComm;
        MPI_Comm_split(mpi->GetWorld(), color, mpi->GetNodeIndex(), &subComm);
        Initialize(subComm);

        return true;
    }
    MPIKernel* MPISub::GetParent(){
        return m_Parent;
    }
    void MPISub::CleanSub(){
        if( IsNodeValid() ){
            MPI_Comm_free(&m_World);
        }
    }
}
