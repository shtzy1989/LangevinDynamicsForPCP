#include "streamer.h"
#include "MOLUTILITY.h"

bool StreamerReader::Open(text filename, std::ios::openmode mode) {
    file.open(filename.c_str(), mode);
    if (file.is_open()) {
        return true;
    }
    else
    {
        return false;
    }
};

void StreamerReader::MoveToBeginning() {
    file.seekg(0, file.beg);
};

void StreamerReader::MoveTo(long long int pos){
    file.seekg(pos);
};

void StreamerReader::MoveToFromCurrent(long long int pos){
    file.seekg(pos, std::ios::cur);
};

long long int StreamerReader::GetFileLength(){
    file.seekg(0, file.end);
    long long int length = file.tellg();
    file.seekg(0);
    return length;
};

long long int StreamerReader::GetPosition(){
    long long int length = file.tellg();
    return length;
};

bool StreamerReader::NextLineAvailable() {
    return file.good();
};

bool StreamerReader::EndOfFile(){
    return file.eof();
};

bool StreamerReader::NextEndOfFile(){
    return file.peek() == file.eof();
};

bool StreamerReader::IsAtTheEnd(){
    long long pos = file.tellg();
    long long length = GetFileLength();
    file.seekg(pos, std::ios::beg);
    return pos == length;
};

double StreamerReader::ReadDouble() {
    double data;
    file.read((char*)&data, sizeof(double));
    return data;
};

float StreamerReader::ReadFloat() {
    float data;
    file.read((char*)&data, sizeof(float));
    return data;
};

int StreamerReader::ReadInt() {
    int data;
    file.read((char*)&data, sizeof(int));
    return data;
};

unsigned int StreamerReader::ReadUInt() {
    unsigned int data;
    file.read((char*)&data, sizeof(unsigned int));
    return data;
};

long StreamerReader::ReadLong() {
    long data;
    file.read((char*)&data, sizeof(long));
    return data;
};

text StreamerReader::ReadText() {
    int length;
    file.read((char*)&length, sizeof(int));
    char* buffer = new char[length + 1];
    file.read(buffer, length);
    buffer[length] = '\0';
    text data = buffer;
    delete[] buffer;
    return data;
}

text StreamerReader::ReadLine() {
    std::string buffer;
    std::getline(file, buffer);
    text data = buffer;
    return data;
}

Real StreamerReader::ReadReal() {
    Real data;
    file.read((char*)&data, sizeof(Real));
    return data;
}

size_t StreamerReader::ReadSizeT() 
{
    size_t data;
    file.read((char*)&data, sizeof(size_t));
    return data;
};

bool StreamerReader::ReadBool()
{
    bool data;
    file.read((char*)&data, sizeof(bool));
    return data;
};

void StreamerReader::ReadChunk(void* dest, size_t size){
    file.read((char*)dest, size);
};

void StreamerReader::Close()
{
    file.close();
};

void StreamerReader::Clear(){
    file.clear();
};

bool StreamerWriter::Open(text filename, std::ios::openmode mode, bool bBackUp) {
    if( bBackUp ){
        MOLUTILITY::findNextBackUpFileName(filename.string());
    }
    file.open(filename.c_str(), mode);
    size = 0;
    if (file.is_open()) {
        return true;
    }
    else
    {
        return false;
    }
};

void StreamerWriter::MoveToBeginning() {
    file.seekp(0);
};

void StreamerWriter::WriteDouble(double data) {
    file.write((char*)&data, sizeof(double));
    size += sizeof(double);
};

void StreamerWriter::WriteFloat(float data) {
    file.write((char*)&data, sizeof(float));
    size += sizeof(float);
};

void StreamerWriter::WriteInt(int data) {
    file.write((char*)&data, sizeof(int));
    size += sizeof(int);
};

void StreamerWriter::WriteUInt(unsigned int data) {
    file.write((char*)&data, sizeof(unsigned int));
    size += sizeof(unsigned int);
};

void StreamerWriter::WriteLong(long data) {
    file.write((char*)&data, sizeof(long));
    size += sizeof(long);
};

void StreamerWriter::WriteText(text data) {
    int length = data.size();
    file.write((char*)&length, sizeof(int));
    file.write(&data[0], length);
    size += sizeof(int) + length;
}

void StreamerWriter::WriteLine(text data) {
    file << data.string() << std::endl;
}

void StreamerWriter::WriteReal(Real data) {
    file.write((char*)&data, sizeof(Real));
    size += sizeof(Real);
}

void StreamerWriter::WriteSizeT(size_t data) {
    file.write((char*)&data, sizeof(size_t));
    size += sizeof(size_t);
}

void StreamerWriter::WriteBool(bool data) {
    file.write((char*)&data, sizeof(bool));
    size += sizeof(bool);
}

void StreamerWriter::WriteChunk(void* dest, size_t size){
    file.write((char*)dest, size);
    this->size += size;
};

void StreamerWriter::Close()
{
    file.close();
};

size_t StreamerWriter::GetCurrentFileLength(){
    return size;
};

bool ReaderTextMatrixFile::Read(text filename){
    StreamerReader reader;
    if( reader.Open(filename) ){
        m_Data.clear();
        while( reader.NextLineAvailable() ){
            text line = reader.ReadLine();
            std::vector<text> word = line.split();
            if( word.size() != 0 ){
                m_Data.push_back(word);
            }
        }
        reader.Close();
        return true;
    }else{
        return false;
    }
};

bool ReaderNumericalTableFile::Read(text filename, text marker, int nheader, int nscale){
    ReaderTextMatrixFile file;
    if( file.Read(filename) ){
        int nword = -1;
        for(int i=0;i<file.m_Data.size();i++){
            if( file.m_Data[i][0][0] != '#' ){
                nword = file.m_Data[i].size();
                break;
            }
        }
        std::vector<std::vector<text> > data;
        for(int i=0;i<file.m_Data.size();i++){
            if( file.m_Data[i].size() == nword ){
                data.push_back(file.m_Data[i]);
            }
        }

        if( nheader == -1 ){
            nheader = 0;
            for(int i=0;i<data.size();i++){
                if( data[i][0][0] == '#' ){
                    nheader++;
                }else{
                    break;
                }
            }
        }
        if( nheader == 0 && nscale == -1 ){
            fprintf(stderr, "Error: can't determine number of columns for scale given no header lines\n");
            return false;
        }
        if( nscale == -1 ){
            nscale = 0;
            for(int i=0;i<data[0].size();i++){
                if( data[0][i][0] == '#' ){
                    nscale++;
                }else{
                    break;
                }
            }
        }

        int ncolumn = data[0].size() - nscale;
        int nrow = data.size() - nheader;

        VectorAuxiliary::AllocateVector2D(&m_Data, nrow, ncolumn);
        VectorAuxiliary::AllocateVector2D(&m_Header, nheader, ncolumn);
        VectorAuxiliary::AllocateVector2D(&m_Scale, nrow, nscale);

        m_NumberOfHeader = nheader;
        m_NumberOfScale = nscale;
        m_NumberOfColumn = ncolumn;
        m_NumberOfRow = nrow;

        // fprintf(stderr, "Numebr of header %10d\n", m_NumberOfHeader);
        // fprintf(stderr, "Numebr of scale  %10d\n", m_NumberOfScale);
        // fprintf(stderr, "Numebr of column %10d\n", m_NumberOfColumn);
        // fprintf(stderr, "Numebr of row    %10d\n", m_NumberOfRow);

        // copy data
        for(int j=0;j<m_NumberOfRow;j++){
            for(int i=0;i<m_NumberOfColumn;i++){
                if( data[j+m_NumberOfHeader][i+m_NumberOfScale].ContentEqualCaseInsensitive("nan") ||
                    data[j+m_NumberOfHeader][i+m_NumberOfScale].ContentEqualCaseInsensitive("-nan") ){
                    m_Data[j][i] = 0.0;
                }else{
                    m_Data[j][i] = data[j+m_NumberOfHeader][i+m_NumberOfScale].to_double();
                }
            }
        }

        // copy header
        for(int j=0;j<m_NumberOfHeader;j++){
            for(int i=0;i<m_NumberOfColumn;i++){
                if( data[j][i+m_NumberOfScale].ContentEqualCaseInsensitive("nan") ||
                    data[j][i+m_NumberOfScale].ContentEqualCaseInsensitive("-nan") ){
                    m_Header[j][i] = 0.0;
                }else{
                    m_Header[j][i] = data[j][i+m_NumberOfScale].to_double();
                }
            }
        }

        // copy scale
        for(int j=0;j<m_NumberOfRow;j++){
            for(int i=0;i<m_NumberOfScale;i++){
                if( data[j+m_NumberOfHeader][i].ContentEqualCaseInsensitive("nan") ||
                    data[j+m_NumberOfHeader][i].ContentEqualCaseInsensitive("-nan") ){
                    m_Scale[j][i] = 0.0;
                }else{
                    m_Scale[j][i] = data[j+m_NumberOfHeader][i].to_double();
                }
            }
        }

        return true;
    }else{
        return false;
    }
};

bool ReaderNumericalTableFile::ReadFixedWith(text filename, text marker, int nheader, int nscale, int scaleWidth, int dataWidth){
    StreamerReader reader;
    if( reader.Open(filename) ){
        std::vector<text> lines;
        while( reader.NextLineAvailable() ){
            text line = reader.ReadLine();
            if( !line.isEmpty() ){
                lines.push_back(line);
            }
        }
        std::vector<std::vector<text> > data;
        if( lines.size() == 0 ){
            fprintf(stderr, "Error: no lines in file\n");
            reader.Close();
            return false;
        }
        int lineLength = lines[0].size();
        int iLine = 0;
        int ndata = (lineLength - nscale * scaleWidth) / dataWidth;
        if( ndata * dataWidth != lineLength - nscale * scaleWidth ){
            fprintf(stderr, "Error: line is not of fixed width\n");
            return false;
        }
        for(int i=0;i<lines.size();i++){
            if( lines[i].size() == lineLength ){
                data.resize(data.size() + 1);
                int pos = 0;
                for(int j=0;j<nscale;j++){
                    data.back().push_back(lines[i].substr(pos, scaleWidth));
                    data.back().back() = data.back().back().replaceAll(marker, " ");
                    pos += scaleWidth;
                }
                for(int j=0;j<ndata;j++){
                    data.back().push_back(lines[i].substr(pos, dataWidth));
                    data.back().back() = data.back().back().replaceAll(marker, " ");
                    pos += dataWidth;
                }
                iLine++;
            }
        }

        int ncolumn = data[0].size() - nscale;
        int nrow = data.size() - nheader;

        VectorAuxiliary::AllocateVector2D(&m_Data, nrow, ncolumn);
        VectorAuxiliary::AllocateVector2D(&m_Header, nheader, ncolumn);
        VectorAuxiliary::AllocateVector2D(&m_Scale, nrow, nscale);

        m_NumberOfHeader = nheader;
        m_NumberOfScale = nscale;
        m_NumberOfColumn = ncolumn;
        m_NumberOfRow = nrow;

        // fprintf(stderr, "Numebr of header %10d\n", m_NumberOfHeader);
        // fprintf(stderr, "Numebr of scale  %10d\n", m_NumberOfScale);
        // fprintf(stderr, "Numebr of column %10d\n", m_NumberOfColumn);
        // fprintf(stderr, "Numebr of row    %10d\n", m_NumberOfRow);

        // copy data
        for(int j=0;j<m_NumberOfRow;j++){
            for(int i=0;i<m_NumberOfColumn;i++){
                if( data[j+m_NumberOfHeader][i+m_NumberOfScale].ContentEqualCaseInsensitive("nan") ||
                    data[j+m_NumberOfHeader][i+m_NumberOfScale].ContentEqualCaseInsensitive("-nan") ){
                    m_Data[j][i] = 0.0;
                }else{
                    m_Data[j][i] = data[j+m_NumberOfHeader][i+m_NumberOfScale].to_double();
                }
            }
        }

        // copy header
        for(int j=0;j<m_NumberOfHeader;j++){
            for(int i=0;i<m_NumberOfColumn;i++){
                if( data[j][i+m_NumberOfScale].ContentEqualCaseInsensitive("nan") ||
                    data[j][i+m_NumberOfScale].ContentEqualCaseInsensitive("-nan") ){
                    m_Header[j][i] = 0.0;
                }else{
                    m_Header[j][i] = data[j][i+m_NumberOfScale].to_double();
                }
            }
        }

        // copy scale
        for(int j=0;j<m_NumberOfRow;j++){
            for(int i=0;i<m_NumberOfScale;i++){
                if( data[j+m_NumberOfHeader][i].ContentEqualCaseInsensitive("nan") ||
                    data[j+m_NumberOfHeader][i].ContentEqualCaseInsensitive("-nan") ){
                    m_Scale[j][i] = 0.0;
                }else{
                    m_Scale[j][i] = data[j+m_NumberOfHeader][i].to_double();
                }
            }
        }

        reader.Close();
        return true;
    }else{
        return false;
    }
};

bool ReaderNumericalTableFile::Write(text filename, bool bTranspose){
    FILE* fout = fopen(filename.c_str(), "w+");
    if( fout == 0 ){
        return false;
    }else{
        if( bTranspose ){
            for(int i=0;i<m_NumberOfScale;i++){
                for(int j=0;j<m_NumberOfHeader;j++){
                    fprintf(fout, "%23s ", "#");
                }
                for(int j=0;j<m_NumberOfRow;j++){
                    fprintf(fout, "%23.16E ", m_Scale[j][i]);
                }
                fprintf(fout, "\n");
            }
            for(int i=0;i<m_NumberOfColumn;i++){
                for(int j=0;j<m_NumberOfHeader;j++){
                    fprintf(fout, "%23.16E ", m_Header[j][i]);
                }
                for(int j=0;j<m_NumberOfRow;j++){
                    fprintf(fout, "%23.16E ", m_Data[j][i]);
                }
                fprintf(fout, "\n");
            }
        }else{
            for(int j=0;j<m_NumberOfHeader;j++){
                for(int i=0;i<m_NumberOfScale;i++){
                    fprintf(fout, "%23s ", "#");
                }
                for(int i=0;i<m_NumberOfColumn;i++){
                    fprintf(fout, "%23.16E ", m_Header[j][i]);
                }
                fprintf(fout, "\n");
            }
            for(int j=0;j<m_NumberOfRow;j++){
                for(int i=0;i<m_NumberOfScale;i++){
                    fprintf(fout, "%23.16E ", m_Scale[j][i]);
                }
                for(int i=0;i<m_NumberOfColumn;i++){
                    fprintf(fout, "%23.16E ", m_Data[j][i]);
                }
                fprintf(fout, "\n");
            }
        }
        fclose(fout);
        return true;
    }
}

void ReaderNumericalTableFile::CopyFrom(ReaderNumericalTableFile* rhs){
    this->m_Data = rhs->m_Data;
    this->m_Header = rhs->m_Header;
    this->m_NumberOfColumn = rhs->m_NumberOfColumn;
    this->m_NumberOfHeader = rhs->m_NumberOfHeader;
    this->m_NumberOfRow = rhs->m_NumberOfRow;
    this->m_NumberOfScale = rhs->m_NumberOfScale;
    this->m_Scale = rhs->m_Scale;
}

bool TextFileWriter::Open(text filename, bool bBackUp){
    if( bBackUp ){
        MOLUTILITY::findNextBackUpFileName(filename.string());
    }
    m_Fout = fopen(filename.c_str(), "w+");
    return m_Fout != 0;
};
void TextFileWriter::Write(const char *format, ...){
    va_list args;
    va_start(args, format);
    vfprintf(m_Fout, format, args);
    va_end(args);
};
void TextFileWriter::WriteFileAndScreen(const char *format, ...){
    va_list args;
    va_start(args, format);
    vfprintf(m_Fout, format, args);
    vfprintf(stderr, format, args);
    va_end(args);
};
void TextFileWriter::Close(){
    fclose(m_Fout);
};
FILE* TextFileWriter::GetFile(){
    return m_Fout;
}