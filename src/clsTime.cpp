#include "clsTime.h"

void clsSystemTime::Print_Time(){
    time_t _time;
    time(&_time);
    fprintf(stderr, "%s\n", ctime(&_time));
};
void clsSystemTime::Print_Time(char* out, const char* format){
    time_t _time;
    time(&_time);
    sprintf(out, format, ctime(&_time));
};
void clsSystemTime::Start(){
    time(&timer);
};
double clsSystemTime::Read(){
    time_t _time;
    time(&_time);
    return difftime(_time, timer);
};
double clsSystemTime::Stop(){
    time_t _time;
    time(&_time);
    double diff = difftime(_time, timer);
    timer = _time;
    return diff;
};
void clsSystemTime::Print_Time_Diff(double second, char* out, const char* format){
    double time_out = second;
    std::string unit = "s";
    double time_in = time_out / 60.0;
    if( time_in > 1.0 ){
        time_out = time_in;
        unit = "m";
        time_in = time_out / 60.0;
        if( time_in > 1.0 ){
            time_out =  time_in;
            unit = "h";
            time_in = time_out / 24.0;
            if( time_in > 1.0 ){
                time_out = time_in;
                unit = "d";
            }
        }
    }
    sprintf(out, format, time_out, unit.c_str());
};

void clsSystemTimeCPP::AddEntry(std::string title){
    m_entryContent.push_back(0);
    m_entryTitle[title] = m_entryContent.size();
};

void clsSystemTimeCPP::StopToEntry(int ID){
    if( ID >= 0 && ID < m_entryContent.size() ){
        m_entryContent[ID] = Stop();
    }
};

void clsSystemTimeCPP::StopToEntry(std::string title){
    StopToEntry(m_entryTitle[title] - 1);
};

std::string clsSystemTimeCPP::ToString(int ID, const char *format){
    if( ID >= 0 && ID < m_entryContent.size() ){
        char buffer[256];
        Print_Time_Diff(m_entryContent[ID], buffer, "%f %s");
        return buffer;
    }else{
        return "";
    }
};

std::string clsSystemTimeCPP::ToString(std::string title, const char *format){
    return ToString(m_entryTitle[title] - 1, format);
};

double clsSystemTimeCPP::GetTime(std::string title){
    return GetTime(m_entryTitle[title] - 1);
};
double clsSystemTimeCPP::GetTime(int ID){
    if( ID >= 0 && ID < m_entryContent.size() ){
        return m_entryContent[ID];
    }else{
        return 0.0;
    }
};

text clsSystemTimeCPP::Print(double value, const char *format){
    char buffer[256];
    Print_Time_Diff(value, buffer, "%f %s");
    return buffer;
}

// ========================================
LoopTimer::LoopTimer(){
    m_FormatValue = "%6.2f";
    m_FormatUnit = "%1s";
    SetFormat();
};

void LoopTimer::Start(){
    m_Timer.Start();
    strcpy(m_TimeToGoString, "Unknown");
    strcpy(m_TimeTotalString, "Unknown");
    sprintf(m_Buffer, "%s/%s", m_TimeToGoString, m_TimeTotalString);
    m_TimeTick = m_Timer.Read();
};

void LoopTimer::Update(int current, int total){
    double timeElapsed = m_Timer.Read();
    double timeToGo = timeElapsed / double(current + 1) * (total - current - 1);
    double timeTotal = timeElapsed / double(current + 1) * (total);
    m_Timer.Print_Time_Diff(timeToGo, m_TimeToGoString, m_Format.c_str());
    m_Timer.Print_Time_Diff(timeTotal, m_TimeTotalString, m_Format.c_str());
    sprintf(m_Buffer, "%s/%s", m_TimeToGoString, m_TimeTotalString);
};

void LoopTimer::Update(long current, long total){
    double timeElapsed = m_Timer.Read();
    double timeToGo = timeElapsed / double(current + 1) * (total - current - 1);
    double timeTotal = timeElapsed / double(current + 1) * (total);
    m_Timer.Print_Time_Diff(timeToGo, m_TimeToGoString, m_Format.c_str());
    m_Timer.Print_Time_Diff(timeTotal, m_TimeTotalString, m_Format.c_str());
    sprintf(m_Buffer, "%s/%s", m_TimeToGoString, m_TimeTotalString);
};

void LoopTimer::Update(long long current, long long total){
    double timeElapsed = m_Timer.Read();
    double timeToGo = timeElapsed / double(current + 1) * (total - current - 1);
    double timeTotal = timeElapsed / double(current + 1) * (total);
    m_Timer.Print_Time_Diff(timeToGo, m_TimeToGoString, m_Format.c_str());
    m_Timer.Print_Time_Diff(timeTotal, m_TimeTotalString, m_Format.c_str());
    sprintf(m_Buffer, "%s/%s", m_TimeToGoString, m_TimeTotalString);
};

void LoopTimer::Update(size_t current, size_t total){
    double timeElapsed = m_Timer.Read();
    double timeToGo = timeElapsed / double(current + 1) * (total - current - 1);
    double timeTotal = timeElapsed / double(current + 1) * (total);
    m_Timer.Print_Time_Diff(timeToGo, m_TimeToGoString, m_Format.c_str());
    m_Timer.Print_Time_Diff(timeTotal, m_TimeTotalString, m_Format.c_str());
    sprintf(m_Buffer, "%s/%s", m_TimeToGoString, m_TimeTotalString);
};

text LoopTimer::GetText(){
    return m_Buffer;
};

char* LoopTimer::GetPChar(){
    return m_Buffer;
};

bool LoopTimer::Tick(double value){
    double current = m_Timer.Read();
    if( current - m_TimeTick >= value ){
        m_TimeTick = current;
        return true;
    }else{
        return false;
    }
}
void LoopTimer::SetValueFormat(text format){
    m_FormatValue = format;
    SetFormat();
}
void LoopTimer::SetUnitFormat(text format){
    m_FormatUnit = format;
    SetFormat();
}
void LoopTimer::SetFormat(){
    m_Format = m_FormatValue + " " + m_FormatUnit;
}
