#ifndef MD_STREAMER_H
#define MD_STREAMER_H

#include <fstream>
#include "text.h"
#include "mathdef.h"

class StreamerReader {
public:
    std::ifstream file;
public:
    bool Open(text filename, std::ios::openmode mode = std::ios::in);
    void MoveToBeginning();
    void MoveTo(long long int pos);
    void MoveToFromCurrent(long long int pos);
    long long int GetFileLength();
    long long int GetPosition();
    bool NextLineAvailable();
    bool EndOfFile();
    bool NextEndOfFile();
    bool IsAtTheEnd();
    double ReadDouble();
    float ReadFloat();
    int ReadInt();
    unsigned int ReadUInt();
    long ReadLong();
    text ReadText();
    text ReadLine();
    Real ReadReal();
    size_t ReadSizeT();
    bool ReadBool();
    template<typename T> std::vector<T> ReadVector(size_t size){
        std::vector<T> data(size);
        file.read((char*)&data[0], sizeof(T) * size);
        return data;
    };
    template<typename T> void ReadVector(T* pointer, size_t size){
        file.read((char*)pointer, sizeof(T) * size);
    };
    template<typename T> T Read(){
        T data;
        file.read((char*)&data, sizeof(T));
        return data;
    };
    void ReadChunk(void* dest, size_t size);
    void Close();
    void Clear();
};

class StreamerWriter {
protected:
    std::ofstream file;
    size_t size;
public:
    bool Open(text filename, std::ios::openmode mode = std::ios::out, bool bBackUp = false);
    void MoveToBeginning();
    void WriteDouble(double data);
    void WriteFloat(float data);
    void WriteInt(int data);
    void WriteUInt(unsigned int data);
    void WriteLong(long data);
    void WriteText(text data);
    void WriteLine(text data);
    void WriteReal(Real data);
    void WriteSizeT(size_t data);
    void WriteBool(bool data);
    void WriteChunk(void* dest, size_t size);
    template<typename T> void WriteVector(std::vector<T>* data){
        file.write((char*)&(*data)[0], sizeof(T) * data->size());
        size += sizeof(T) * data->size();
    };
    template<typename T> void Write(T data){
        file.write((char*)&data, sizeof(T));
        size += sizeof(T);
    };
    void Close();
    size_t GetCurrentFileLength();
    static std::string GetReadableFileLength(size_t length){
        std::string result = "0";
        char sunit[5] = { 'T', 'G', 'M', 'K', 'B'};
        double value = 0.0;
        size_t unit[6];
        unit[5] = 0;
        unit[4] = length % 1024; length /= 1024;
        unit[3] = length % 1024; length /= 1024;
        unit[2] = length % 1024; length /= 1024;
        unit[1] = length % 1024; length /= 1024;
        unit[0] = length;
        for(int i=0;i<5;i++){
            if( unit[i] ){
                value = double(unit[i]) + double(unit[i+1]) / 1024.0;
                char buffer[256];
                sprintf(buffer, "%6.2f %c", value, sunit[i]);
                result = buffer;
                break;
            }
        }
        return result;
    };
};

class ReaderTextMatrixFile{
public:
    std::vector<std::vector<text> > m_Data;
public:
    bool Read(text filename);
};

class ReaderNumericalTableFile{ // Row major
public:
    std::vector<std::vector<double> > m_Header;
    std::vector<std::vector<double> > m_Scale;
    std::vector<std::vector<double> > m_Data;
    int m_NumberOfHeader;
    int m_NumberOfScale;
    int m_NumberOfColumn;
    int m_NumberOfRow;
public:
    bool ReadFixedWith(text filename, text marker, int nheader, int nscale, int scaleWidth = 14, int dataWidth = 24);
    bool Read(text filename, text marker = "#", int nheader = -1, int nscale = -1);
    bool Write(text filename, bool bTranspose = false);
    void CopyFrom(ReaderNumericalTableFile* rhs);
};

class TextFileWriter{
protected:
    FILE* m_Fout;
public:
    bool Open(text filename, bool bBackUp = false);
    void Write(const char *format, ...);
    void WriteFileAndScreen(const char *format, ...);
    void Close();
    FILE* GetFile();
};


#endif

