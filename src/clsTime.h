#ifndef MD_TIME_H
#define MD_TIME_H

#include "cppheader.h"
#include "text.h"

class clsSystemTime{
private:
    time_t timer;
    struct tm t2k;
public:
    void Print_Time();
    void Print_Time(char* out, const char* format);
    void Start();
    double Read();
    double Stop();
    static void Print_Time_Diff(double second, char* out, const char* format);
};

class clsSystemTimeCPP : public clsSystemTime{
private:
    std::vector<double> m_entryContent;
    std::map<std::string, int> m_entryTitle;
public:
    void AddEntry(std::string title);
    void StopToEntry(int ID);
    void StopToEntry(std::string title);
    std::string ToString(int ID, const char *format = "%f %s");
    std::string ToString(std::string title, const char *format = "%f %s");
    double GetTime(std::string title);
    double GetTime(int ID);
    static text Print(double value, const char *format = "%f %s");
};

class LoopTimer{
private:
    clsSystemTimeCPP m_Timer;
    char m_TimeToGoString[256];
    char m_TimeTotalString[256];
    char m_Buffer[256];
    double m_TimeTick;
    text m_FormatValue;
    text m_FormatUnit;
    text m_Format;
public:
    LoopTimer();
public:
    void Start();
    void Update(int current, int total);
    void Update(long current, long total);
    void Update(long long current, long long total);
    void Update(size_t current, size_t total);
    text GetText();
    char* GetPChar();
    bool Tick(double value);
    void SetValueFormat(text format);
    void SetUnitFormat(text format);
protected:
    void SetFormat();
};

#endif
