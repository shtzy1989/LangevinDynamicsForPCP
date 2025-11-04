#ifndef MD_UTILITY_STRING_H
#define MD_UTILITY_STRING_H

#include "cppheader.h"

class text{
private:
    std::string data;
public:
#ifdef FMT
    template<typename... Args>
    text(const std::string& format, Args&&... args){
        data = fmt::format(format, std::forward<Args>(args)...);
    };
#endif
    text() {};
    text(const char *p);
    text(char *p);
    text(std::string p);
    text(bool p);
    text(int p);
    text(uint p);
    text(size_t p);
    text(float p);
    text(long p);
    text(long long p);
    text(double p);
    text(char p);
    text(char* p, int len);
    text(const char *format, int p);
    text(const char *format, uint p);
    text(const char *format, size_t p);
    text(const char *format, float p);
    text(const char *format, double p);
    text(const char *format, char *p);
    ~text() {};

    text& operator = (const char *p);
    text& operator = (char *p);
    text& operator = (const std::string& p);
    text& operator = (bool p);
    text& operator = (const text& rhs);
    char& operator[](size_t p);
    const char& operator[](size_t p) const;
    text& operator += (const text& rhs);

    void push_back(char c);
    void insert(size_t p, char c);
    void insert(size_t p, std::string c);
    void erase(size_t p);
    void erase(size_t p, size_t len);
    size_t size();
    void resize(size_t size);

    text substr(size_t p, size_t stop = std::string::npos) const;
    text slice(int start, int end = INT_MAX, int stride = 1) const;
    text left(size_t len) const;
    text right(size_t len) const;
    bool leftEqual(text value) const;
    bool rightEqual(text value) const;    
    double to_double() const;
    int to_int() const;
    size_t to_unsigned_long() const;
    float to_float() const;
    bool to_boolean() const;
    long to_long() const;
    long long to_longlong() const;
    int letterToInt() const;

    std::vector<text> split(text mask = " \t");
    std::vector<int> splitToInt(text mask = " \t");
    std::vector<float> splitToFloat(text mask = " \t");
    std::vector<double> splitToDouble(text mask = " \t");
    std::vector<text> splitAdvanced(Array<text> leftBracket, Array<text> rightBracket, text mask = " \t");
    std::vector<text> splitAdvanced(text leftBracket, text rightBracket, text mask = " \t");
    text parsePath();
    text parseFileName();
    text parseRemoveComment(text symbol);
    text parseGetComment(text symbol);
    text parseExtension();
    text removeExtension();
    std::vector<text> parseBracket(text left, text right); // can find one multi-character bracket
    std::vector<text> parseBracket(text symbol); // can find one multi-character bracket
    text parseFirstBracket(text left, text right); // can find one multi-character bracket
    text parseFirstBracket(text symbol); // can find one multi-character bracket
    std::vector<text> parseOuterBracket(text left, text right); // used to find bracket that has inner brackets
    text replaceOuterBracket(text left, text right, std::vector<text> value); // used to find bracket that has inner brackets
    text replaceAll(text word, text newword);
    text trim(size_t len);
    text removeSpace();
    text removeSpaceBeforeAfter();
    bool isEmpty() const;
    bool isDouble() const;
    bool isFloat() const;
    bool isInt() const;
    bool isLong() const;
    bool isLongLong() const;
    bool isBoolean() const;
    std::vector<size_t> Find(text word);
    bool Contain(text word) const;
    bool ContentEqual(text rhs, text mask = " \t");
    bool ContentEqualCaseInsensitive(text rhs, text mask = " \t");
    bool CompareWithAsteriskMask(text mask);
    text GetEqualSignValue(text lhs, text equalSign = "=", text delimiter = ", \t");
    std::vector<std::pair<text, text> > GetEqualSignValueAll(text equalSign = "=", text delimiter = ", \t");
    int numberOfAppearances(text word) const;

    size_t find_first_of(text symbol) { return data.find_first_of(symbol.data); };
    size_t find_first_not_of(text symbol) { return data.find_first_not_of(symbol.data); };
    size_t find_last_of(text symbol) { return data.find_last_of(symbol.data); };
    size_t find_last_not_of(text symbol) { return data.find_last_not_of(symbol.data); };
    size_t find(char c, size_t pos){ return data.find(c, pos); };
    size_t find(const char* c, size_t pos){ return data.find(c, pos); };
    size_t find(const std::string& c, size_t pos){ return data.find(c, pos); };
    size_t find(const char* c, size_t pos, size_t n){ return data.find(c, pos, n); };

    const char *c_str() const;
    const std::string string() const;

    static text FromTextArray(std::vector<text> *word, int first = 0, int last = -1);
    static text FromDoubleArray(std::vector<double> *word, text format = "%13.7f", int first = 0, int last = -1);
    static text FromIntArray(std::vector<int> *word, text format = "%13d", int first = 0, int last = -1);
    static text GetLine(std::ifstream *fin);
    static text FromByte2Readable(double value);
    static text FromMilliSecond2Readable(double value);
    static text FromSecond2Readable(double value);

    static text formatC(const char *format, ...);

    friend std::ostream& operator<<(std::ostream& lhs, const text& rhs);
}; 

text operator + (const text& lhs, const text& rhs);
bool operator < (const text& lhs, const text& rhs);
bool operator > (const text& lhs, const text& rhs);
bool operator == (const text& lhs, const text& rhs);
bool operator != (const text& lhs, const text& rhs);
text operator * (const text& lhs, const int& rhs);
std::ostream& operator<<(std::ostream& lhs, const text& rhs);

#endif
