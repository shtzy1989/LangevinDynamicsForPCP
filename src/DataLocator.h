#ifndef MD_DATA_LOCATOR_H
#define MD_DATA_LOCATOR_H
#include "cppheader.h"
#include "streamer.h"

namespace dataLocator{
    template<typename T, typename T2>
    class Element2Tag{
    protected:
        T m_Value[2];
        std::vector<T2> m_Tag;
    public:
        Element2Tag(T value0, T value1, std::vector<T2> tag, bool bOrder){
            if( bOrder ){
                if( value0 < value1 ){
                    m_Value[0] = value0;
                    m_Value[1] = value1;
                }else{
                    m_Value[0] = value1;
                    m_Value[1] = value0;
                }
            }else{
                m_Value[0] = value0;
                m_Value[1] = value1;
            }
            m_Tag = tag;
        }
        virtual T Value(int index) const{
            return m_Value[index];
        }
        virtual T2 Tag(int index) const{
            return m_Tag[index];
        }
        size_t TagSize(){
            return m_Tag.size();
        }
        virtual bool operator<(const Element2Tag& rhs) const{
            return ( m_Value[0] == rhs.m_Value[0] && 
                      m_Value[1] < rhs.m_Value[1] ) || 
                      m_Value[0] < rhs.m_Value[0];
        }
    };

    template<typename T, typename T2>
    class Element2Tag1 : public Element2Tag<T, T2>{
    public:
        Element2Tag1(T value0, T value1, T2 tag, bool bOrder) : Element2Tag<T, T2>(value0, value1, {tag}, bOrder){
        };
    };

    template<typename T, typename T2>
    class Element2Tag1Order : public Element2Tag<T, T2>{
    public:
        Element2Tag1Order(T value0, T value1, T2 tag) : Element2Tag<T, T2>(value0, value1, {tag}, true){
        };
    };

    template<typename T, typename T2>
    class Element2Tag1Directional : public Element2Tag<T, T2>{
    public:
        Element2Tag1Directional(T value0, T value1, T2 tag) : Element2Tag<T, T2>(value0, value1, {tag}, false){
        };
    };

    template<typename T>
    class Element2 : public Element2Tag<T, int>{
    public:
        Element2(T value0, T value1, bool bOrder) : Element2Tag<T, int>(value0, value1, {}, bOrder){
        };
    };

    template<typename T>
    class Element2Order : public Element2Tag<T, int>{
    public:
        Element2Order(T value0, T value1) : Element2Tag<T, int>(value0, value1, {}, true){
        };
    };

    template<typename T>
    class Element2Directional : public Element2Tag<T, int>{
    public:
        Element2Directional(T value0, T value1) : Element2Tag<T, int>(value0, value1, {}, false){
        };
    };

    // ============================================================================
    // FiniteInifiniteIntegerTickArray1DElement
    // ============================================================================
    class FiniteInifiniteIntegerTickArray1DElement{
    public:
        virtual void Initialize(int ID, void* userData);
        virtual std::string ToText();
        virtual void Write(StreamerWriter *writer);
        virtual void Read(StreamerReader *reader);
    };

    // ============================================================================
    // FiniteInifiniteIntegerTickArray1D
    // ============================================================================
    template<typename T>
    class FiniteInifiniteIntegerTickArray1D{
    protected:
        std::vector<FiniteInifiniteIntegerTickArray1DElement*> m_Data;
        int m_MinID;
        void* m_UserData;
    public:
        FiniteInifiniteIntegerTickArray1D(){
            m_MinID = 0;
            m_UserData = 0;
        }
        ~FiniteInifiniteIntegerTickArray1D(){
            Release();
        }
    public:
        FiniteInifiniteIntegerTickArray1DElement* operator[](int ID){
            return GetElementByID(ID);
        }
        FiniteInifiniteIntegerTickArray1DElement* GetElementByIndex(int index){
            return m_Data[index];
        }
        FiniteInifiniteIntegerTickArray1DElement* GetElementByID(int ID){
            if( m_Data.size() == 0 ){
                m_MinID = ID;
                m_Data.resize(1, 0);
                m_Data[0] = new T;
                m_Data[0]->Initialize(ID, m_UserData);
            }else if( ID < m_MinID ){
                int sizeNew = m_MinID - ID; 
                m_Data.insert(m_Data.begin(), sizeNew, 0);

                for(int i=ID;i<m_MinID;i++){
                    m_Data[i - ID] = new T;
                    m_Data[i - ID]->Initialize(i, m_UserData);
                }

                m_MinID = ID;
            }else if( ID - m_MinID + 1 > m_Data.size() ){
                int sizeNew = ID - m_MinID + 1 - m_Data.size();
                int oldSize = m_Data.size() + m_MinID;
                m_Data.insert(m_Data.end(), sizeNew, 0);

                for(int i=oldSize;i<=ID;i++){
                    m_Data[i - m_MinID] = new T;
                    m_Data[i - m_MinID]->Initialize(i, m_UserData);
                }
            }
            return m_Data[ID - m_MinID];
        }
        size_t GetSize(){
            return m_Data.size();
        }
        int GetMinID(){
            return m_MinID;
        }
        void Print(FILE* fout = stderr){
            for(int i=0;i<m_Data.size();i++){
                fprintf(fout, "%10d %13s\n", i + m_MinID, m_Data[i]->ToText().c_str());
            }
        }
        void SetUserData(void* userData){
            m_UserData = userData;
        }
        void Write(StreamerWriter *writer){
            writer->WriteSizeT(m_Data.size());
            for(int i=0;i<m_Data.size();i++){
                m_Data[i]->Write(writer);
            }
            writer->WriteInt(m_MinID);
        }
        void Read(StreamerReader *reader){
            Release();
            size_t size = reader->ReadSizeT();
            m_Data.resize(size);
            for(int i=0;i<m_Data.size();i++){
                m_Data[i] = new T;
                m_Data[i]->Read(reader);
            }
            m_MinID = reader->ReadInt();
        }
        void Release(){
            for(int i=0;i<m_Data.size();i++){
                SAFE_DELETE(m_Data[i]);
            }
            m_Data.resize(0);
            m_Data.shrink_to_fit();
            m_MinID = 0;
        }
    };

};

#endif
