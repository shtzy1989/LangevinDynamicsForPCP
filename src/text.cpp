#include "text.h"

text::text(const char *p){
    data = p;
};

text::text(char *p){
    data = p;
};

text::text(std::string p){
    data = p;
};

text::text(bool p){
    data = p ? "true" : "false";
};

text::text(int p){
    data = std::to_string(p);
};

text::text(uint p){
    data = std::to_string(p);
};

text::text(size_t p){
    data = std::to_string(p);
};

text::text(float p){
    data = std::to_string(p);
};

text::text(long p){
    data = std::to_string(p);
};

text::text(long long p){
    data = std::to_string(p);
};

text::text(double p){
    data = std::to_string(p);
};

text::text(char p){
    data = p;
};

text::text(const char *format, int p){
    char *buffer = new char[65536]; 
    sprintf(buffer, format, p);
    data = buffer;
    delete [] buffer;
};

text::text(const char *format, uint p){
    char *buffer = new char[65536];
    sprintf(buffer, format, p);
    data = buffer;
    delete [] buffer;
};

text::text(const char *format, size_t p){
    char *buffer = new char[65536];
    sprintf(buffer, format, p);
    data = buffer;
    delete [] buffer;
};

text::text(const char *format, float p){
    char *buffer = new char[65536];
    sprintf(buffer, format, p);
    data = buffer;
    delete [] buffer;
};

text::text(const char *format, double p){
    char *buffer = new char[65536];
    sprintf(buffer, format, p);
    data = buffer;
    delete [] buffer;
};

text::text(const char *format, char *p){
    char *buffer = new char[65536];
    sprintf(buffer, format, p);
    data = buffer;
    delete [] buffer;
};

text::text(char* p, int len) {
    data = std::string(p, len);
}

text& text::operator = (const char *p){
    data = p;
    return *this;
};

text& text::operator = (char *p){
    data = p;
    return *this;
};

text& text::operator = (const std::string& p){
    data = p;
    return *this;
};

text& text::operator = (bool p){
    data = p ? "true" : "false";
    return *this;
};

text& text::operator = (const text& rhs){
    data = rhs.data;
    return *this;
};

char& text::operator[](size_t p){
    return data[p];
};

const char& text::operator[](size_t p) const{
    return data[p];
};

text& text::operator += (const text& rhs){
    data = data + rhs.data;
    return *this;
};

void text::push_back(char c){
    data.push_back(c);
};

void text::insert(size_t p, char c){
    data.insert(data.begin() + p, c);
};

void text::insert(size_t p, std::string c){
    data = data.substr(0, p) + c + data.substr(p, data.size()-p);
};

void text::erase(size_t p){
    data.erase(p);
};

void text::erase(size_t p, size_t len){
    data.erase(p, len);
};

size_t text::size(){
    return data.size();
};

void text::resize(size_t size){
    data.resize(size);
};

text text::substr(size_t p, size_t len) const{
    return data.substr(p, len);
};

text text::slice(int start, int stop, int stride) const{
    // Adjust negative indices
    int len = data.length();

    if( stride > 0 ){
        if( start < 0 ) start = len + start;
        if( stop < 0 ) stop = len + stop;
        start = std::max(0, start);
        stop = std::min(len, stop);
    }else if( stride == 0 ){
        fprintf(stderr, "Error: stride cannot be 0\n");
        return {};
    }else{
        if( start >= len ) start = len - 1;
        stop = stop < -1 ? -1 : stop; // Allowing stop to be -1 for negative stride
    }

    std::string result;
    if (stride > 0) {
        for (int i = start; i < stop; i += stride) {
            result += data[i];
        }
    } else { // Handle negative stride
        for (int i = start; i > stop; i += stride) {
            result += data[i];
        }
    }
    return result;
}

text text::left(size_t len) const{
    len = std::min<size_t>(len, data.size());
    return data.substr(0, len);
};

text text::right(size_t len) const{
    len = std::min<size_t>(len, data.size());
    return data.substr(std::max<int>(data.size() - len, 0), len);
};

bool text::leftEqual(text value) const{
    return left(value.size()).ContentEqual(value);
};

bool text::rightEqual(text value) const{
    return right(value.size()).ContentEqual(value);
};

double text::to_double() const{
    try{
        return std::stod(data);
    }catch(std::invalid_argument &e){
        fprintf(stderr, "Error: Invalid string or overflow %s\n", data.c_str());
        return 0;
    }
};

int text::to_int() const{
    try{
        return std::stoi(data);
    }catch(std::invalid_argument &e){
        fprintf(stderr, "Error: Invalid string or overflow %s\n", data.c_str());
        return 0;
    }
};

size_t text::to_unsigned_long() const{
    try{
        return std::stoul(data);
    }catch(std::invalid_argument &e){
        fprintf(stderr, "Error: Invalid string or overflow %s\n", data.c_str());
        return 0;
    }
};

long text::to_long() const{
    try{
        return std::stol(data);
    }catch(std::invalid_argument &e){
        fprintf(stderr, "Error: Invalid string or overflow %s\n", data.c_str());
        return 0;
    }
};

long long text::to_longlong() const{
    try{
        return std::stoll(data);
    }catch(std::invalid_argument &e){
        fprintf(stderr, "Error: Invalid string or overflow %s\n", data.c_str());
        return 0;
    }
};

int text::letterToInt() const{
    if( data.size() == 1 && isalpha(data[0]) ){
        if( data[0] >= 65 && data[0] < 91 ){
            return data[0] - 65;
        }else if( data[0] >= 97 && data[0] < 123 ){
            return data[0] - 97;
        }else{
            return 0;
        }
    }else{
        return to_int();
    }
}

float text::to_float() const{
    return std::stof(data);
};

bool text::to_boolean() const{
    return data == "TRUE" || 
           data == "True" || 
           data == "true" || 
           ( isInt() && to_int() != 0 ) ||
           ( isDouble() && to_double() != 0 ) ||
           ( isFloat() && to_float() != 0 );
};

std::vector<text> text::split(text mask){
    std::vector<text> out;
    size_t s = 0;
    for(size_t i=0;i<data.size();i++){
        bool bin = false;
        for(size_t l=0;l<mask.size();l++){
            if( data[i] == mask[l] ){
                bin = true;
                break;
            }
        }
        if( bin ){
            if( i != s ){
                out.push_back(data.substr(s, i - s));
            }
            s = i + 1;
        }else if( i == data.size() - 1 && !bin ){
            out.push_back(data.substr(s, i - s + 1));
        }
    }

    return out;
};

std::vector<int> text::splitToInt(text mask){
    std::vector<text> word = split(mask);
    std::vector<int> result;
    result.reserve(word.size());
    for(int i=0;i<word.size();i++){
        result.push_back(word[i].to_int());
    }
    return result;
};

std::vector<float> text::splitToFloat(text mask){
    std::vector<text> word = split(mask);
    std::vector<float> result;
    result.reserve(word.size());
    for(int i=0;i<word.size();i++){
        result.push_back(word[i].to_float());
    }
    return result;
};

std::vector<double> text::splitToDouble(text mask){
    std::vector<text> word = split(mask);
    std::vector<double> result;
    result.reserve(word.size());
    for(int i=0;i<word.size();i++){
        result.push_back(word[i].to_double());
    }
    return result;
};

std::vector<text> text::splitAdvanced(text leftBracket, text rightBracket, text mask){
    std::vector<text> left, right;
    left.push_back(leftBracket);
    right.push_back(rightBracket);
    return splitAdvanced(left, right, mask);
};

std::vector<text> text::splitAdvanced(Array<text> leftBracket, Array<text> rightBracket, text mask){
    text buffer = data;
    int maxBracket = std::min<int>(leftBracket.size(), rightBracket.size());
    std::vector<text> result;
    size_t s = 0;
    int lbn = -1;
    // find bracket position
    Array<int> bracket(data.size(), -1);
    for(int i=0;i<maxBracket;i++){
        std::vector<size_t> posLeft = Find(leftBracket[i]);
        std::vector<size_t> posRight = Find(rightBracket[i]);
        for(int j=0;j<posLeft.size();j++){
            bracket[posLeft[j]] = i; // left bracket begin
        }
        for(int j=0;j<posRight.size();j++){ 
            bracket[posRight[j] + rightBracket[i].size() - 1] = maxBracket + i; // right bracket end
        }
    }

    for(size_t i=0;i<data.size();i++){
        if( 0 <= bracket[i] && bracket[i] < maxBracket ){ // left bracket
            if( i != s ){
                result.push_back(data.substr(s, i - s));
            }
            s = i;
            lbn = bracket[i];
        }else if( maxBracket <= bracket[i] && bracket[i] < maxBracket * 2 ){
            if( lbn != bracket[i] - maxBracket ){
                // error, bracket does not match
                result.resize(0);
                return result;
            }
            result.push_back(data.substr(s, i - s + 1));
            s = i + 1;
            lbn = -1;
        }else if( lbn == -1 ){
            bool bin = false;
            for(size_t l=0;l<mask.size();l++){
                if( data[i] == mask[l] ){
                    bin = true;
                    break;
                }
            }
            if( bin ){
                if( i != s ){
                    result.push_back(data.substr(s, i - s));
                }
                s = i + 1;
            }else if( i == data.size() - 1 && !bin ){
                result.push_back(data.substr(s, i - s + 1));
            }
        }
    }
    return result;
};

text text::parsePath(){
    size_t pos = data.find_last_of("/");
    if( pos == std::string::npos ) return "";
    else return data.substr(0, pos+1);
};

text text::parseFileName(){
    size_t pos = data.find_last_of("/");
    if( pos == std::string::npos ) return data;
    else return data.substr(pos+1, data.size()-pos-1);
};

text text::parseRemoveComment(text symbol){
    size_t pos = data.find_first_of(symbol.c_str());
    if( pos == std::string::npos ) return data;
    else return data.substr(0, pos);
};

text text::parseGetComment(text symbol){
    size_t pos = data.find_first_of(symbol.c_str());
    if( pos == -1 ) return "";
    else return data.substr(pos+1, data.size() - pos-1);
};

text text::parseExtension(){
    size_t pos = -1;
    for(int i=data.size()-1;i>=0;i--){
        if( data[i] == '.' ){
            pos = i;
            break;
        }
    }
    return data.substr(pos + 1, data.size() - pos - 1);
};

text text::removeExtension(){
    size_t pos = -1;
    for(int i=data.size()-1;i>=0;i--){
        if( data[i] == '.' ){
            pos = i;
            break;
        }
    }
    return data.substr(0, pos);
};

std::vector<text> text::parseBracket(text left, text right){
    std::vector<size_t> posLeft = Find(left);
    std::vector<size_t> posRight = Find(right);
    std::vector<text> result;
    size_t numberOfPair = std::min<size_t>(posLeft.size(), posRight.size());
    for(size_t i=0;i<numberOfPair;i++){
        if( posRight[i] > posLeft[i] )
            result.push_back(data.substr(posLeft[i]+left.size(), posRight[i] - posLeft[i] - left.size()));
    }
    return result;
};

std::vector<text> text::parseBracket(text symbol){
    std::vector<size_t> pos = Find(symbol);
    std::vector<text> result;
    for(size_t i=0;i<pos.size()/2;i++){
        result.push_back(data.substr(pos[i*2]+symbol.size(), pos[i*2+1] - pos[i*2] - symbol.size()));
    }
    return result;
};

text text::parseFirstBracket(text left, text right){
    std::vector<size_t> posLeft = Find(left);
    std::vector<size_t> posRight = Find(right);
    if( posLeft.size() >= 1 && posRight.size() >= 1 ){
        return data.substr(posLeft[0]+left.size(), posRight[0] - posLeft[0] - left.size());
    }else{
        return "";
    }
};

text text::parseFirstBracket(text symbol){
    std::vector<size_t> pos = Find(symbol);
    if( pos.size() >= 2 ){
        return data.substr(pos[0]+symbol.size(), pos[1] - pos[0] - symbol.size());
    }else{
        return "";
    }
};

std::vector<text> text::parseOuterBracket(text left, text right){
    std::vector<std::pair<size_t, size_t> > pairBracket;
    int currentLeft = -1;
    int currentLevel = 0;
    for(size_t i=0;i<data.size();i++){
        if( i + left.size() <= data.size() && data.substr(i, left.size()) == left ){
            // this is a left bracket
            if( currentLevel == 0 ){
                // only capture the first level
                currentLeft = i;
            }
            currentLevel++;
            i += left.size() - 1;
        }else if( i + right.size() <= data.size() && data.substr(i, right.size()) == right ){
            // this is a right bracket
            currentLevel--;
            if( currentLevel == 0 ){
                // this is a level right of the pair
                pairBracket.push_back(std::pair<int, int>(currentLeft, i));
            }
            i += right.size() - 1;
        }
    }
    std::vector<text> result;
    for(int i=0;i<pairBracket.size();i++){
        result.push_back(data.substr(
            pairBracket[i].first + left.size(), 
            pairBracket[i].second - pairBracket[i].first - left.size()));
    }
    return result;
};

text text::replaceOuterBracket(text left, text right, std::vector<text> value){
    std::vector<std::pair<size_t, size_t> > pairBracket;
    int currentLeft = -1;
    int currentLevel = 0;
    for(size_t i=0;i<data.size();i++){
        if( i + left.size() <= data.size() && data.substr(i, left.size()) == left ){
            // this is a left bracket
            if( currentLevel == 0 ){
                // only capture the first level
                currentLeft = i;
            }
            currentLevel++;
            i += left.size() - 1;
        }else if( i + right.size() <= data.size() && data.substr(i, right.size()) == right ){
            // this is a right bracket
            currentLevel--;
            if( currentLevel == 0 ){
                // this is a level right of the pair
                pairBracket.push_back(std::pair<int, int>(currentLeft, i));
            }
            i += right.size() - 1;
        }
    }
    text result = data;
    for(int i=pairBracket.size()-1;i>=0;i--){
        result.data.replace(pairBracket[i].first,
            pairBracket[i].second - pairBracket[i].first + right.size(),
            value[i].c_str());
    }
    return result;
};

text text::replaceAll(text word, text newword){
    text result = data;
    std::vector<size_t> pos = Find(word);
    for(int i=pos.size()-1;i>=0;i--){
        result.data.replace(pos[i], word.size(), newword.c_str());
    }
    return result;
};

text text::trim(size_t len){
    text result = data;
    result.data.resize(std::min<int>(result.data.size(), len));
    return result;
};

text text::removeSpace(){
    text result;
    for(int i=0;i<data.size();i++){
        if( !isspace(data[i]) ){
            result.push_back(data[i]);
        }	
    }
    return result;
};

text text::removeSpaceBeforeAfter(){
    size_t start = -1;
    size_t end = data.size();
    for(size_t i=0;i<data.size();i++){
        if( start == -1 && data[i] != ' ' ){
            start = i;
        }
        if( start != -1 && data[i] == ' ' ){
            end = i;
        }
    }
    if( start == -1 ) return "";
    else return data.substr(start, end - start);
}

bool text::isEmpty() const{
    return data.find_first_not_of(" \t") == std::string::npos;
};

bool text::isDouble() const{
    bool bIsNumber = true;
    try{
        double test = std::stod(data);
        bIsNumber = true;
    }catch( std::exception e ){
        bIsNumber = false;
    }
    return bIsNumber;
};

bool text::isFloat() const{
    bool bIsNumber = true;
    try{
        float test = std::stof(data);
        bIsNumber = true;
    }catch( std::exception e ){
        bIsNumber = false;
    }
    return bIsNumber;
};

bool text::isInt() const{
    bool bIsNumber = true;
    try{
        int test = std::stoi(data);
        bIsNumber = true;
    }catch( std::exception e ){
        bIsNumber = false;
    }
    return bIsNumber;
};

bool text::isBoolean() const{
    if( data == "TRUE" ||
        data == "true" ||
        data == "True" ||
        data == "FALSE" ||
        data == "false" ||
        data == "False" ){
        return true;
    }else{
        return false;
    }
};

bool text::isLong() const{
    bool bIsNumber = true;
    try{
        long test = std::stol(data);
        bIsNumber = true;
    }catch( std::exception e ){
        bIsNumber = false;
    }
    return bIsNumber;
};

bool text::isLongLong() const{
    bool bIsNumber = true;
    try{
        long test = std::stoll(data);
        bIsNumber = true;
    }catch( std::exception e ){
        bIsNumber = false;
    }
    return bIsNumber;
};

std::vector<size_t> text::Find(text word){
    std::vector<size_t> result;
    size_t found = data.find(word.data);
    while( found != std::string::npos ){
        result.push_back(found);
        found = data.find(word.data, found + 1);
    }
    return result;
};

bool text::Contain(text word) const{
    return data.find(word.data) != std::string::npos;
};

bool text::ContentEqual(text rhs, text mask){
    bool bEqual = true;
    std::vector<text> splitLhs = split(mask);
    std::vector<text> splitRhs = rhs.split(mask);
    if( splitLhs.size() != splitRhs.size() ){
        bEqual = false;
    }else{
        for(size_t i=0;i<splitLhs.size();i++){
            if( splitLhs[i] != splitRhs[i] ){
                bEqual = false;
            }
        }
    }
    return bEqual;
};

bool text::ContentEqualCaseInsensitive(text rhs, text mask){
    bool bEqual = true;
    std::vector<text> splitLhs = split(mask);
    std::vector<text> splitRhs = rhs.split(mask);
    if( splitLhs.size() != splitRhs.size() ){
        bEqual = false;
    }else{
        for(size_t i=0;i<splitLhs.size();i++){
            if( splitLhs[i].size() != splitRhs[i].size() ){
                bEqual = false;
            }else{
                for(int j=0;j<splitLhs[i].size();j++){
                    if( tolower(splitLhs[i][j]) != tolower(splitRhs[i][j]) ){
                        bEqual = false;
                    }
                }
            }
        }
    }
    return bEqual;
};

bool text::CompareWithAsteriskMask(text mask){
    std::vector<text> word = mask.split("*");
    size_t searchPoint = 0;
    bool leftAsterisk = mask[0] == '*';
    bool rightAsterisk = mask[mask.size() - 1] == '*';
    size_t pos;
    for(int i=0;i<word.size();i++){
        if( i == 0 &&  !leftAsterisk ){
            pos = data.find(word[i].string(), searchPoint);
            if( pos != 0 ){
                return false;
            }
        }else if( i == word.size() - 1 && !rightAsterisk ){
            pos = 0;
            if( data.size() < word[i].size() || this->right(word[i].size()) != word[i] ){
                return false;
            }
        }else{
            pos = data.find(word[i].string(), searchPoint);
            if( pos == -1 ){
                return false;
            }
        }
        searchPoint = pos + word[i].size();
    }
    return true;
};
text text::GetEqualSignValue(text lhs, text equalSign, text delimiter){
    // auto pos = data.find(lhs.string());
    // if( pos != -1 ){
    //     auto pos2 = -1;
    //     for(int i=pos + lhs.size();i<data.size();i++){
    //         if( data.substr(i, equalSign.size()) == equalSign ){
    //             pos2 = i + equalSign.size();
    //             break;
    //         }
    //     }
    //     if( pos2 != -1 ){
    //         text rhs = data.substr(pos2);
    //         auto word = rhs.split(delimiter);
    //         if( word.size() ){
    //             return word[0];
    //         }else{
    //             return "";
    //         }
    //     }else{
    //         return "";
    //     }
    // }else{
    //     return "";
    // }

    auto word0 = this->split(delimiter);
    std::vector<text> word;
    for(int i=0;i<word0.size();i++){
        if( word0[i].Contain(equalSign) ){
            if( word0[i].left(1).ContentEqualCaseInsensitive(equalSign) ){
                // bugged
            }else if( word0[i].right(1).ContentEqualCaseInsensitive(equalSign) ){
                if( i + 1 < word0.size() && !word0[i+1].Contain(equalSign) ){
                    word.push_back(word0[i] + word0[i+1]);
                    i += 1;
                }else{
                    // bugged
                }
            }else{
                word.push_back(word0[i]);
            }
        }else if( i + 2 < word0.size() &&
            word0[i+1].ContentEqualCaseInsensitive(equalSign) &&
            !word0[i+2].Contain(equalSign) ){
            word.push_back(word0[i] + word0[i+1] + word0[i+2]);
            i += 2;
        }else if( i + 1 < word0.size() && 
            word0[i+1].Contain(equalSign) ){
            word.push_back(word0[i] + word0[i+1]);
            i += 1;
        }
    }
    for(int i=0;i<word.size();i++){
        auto word1 = word[i].split(equalSign);
        if( word1.size() == 2 && word1[0] == lhs ){
            return word1[1];
        }
    }
    return "";
}
std::vector<std::pair<text, text> >  text::GetEqualSignValueAll(text equalSign, text delimiter){
    std::vector<std::pair<text, text> >  result;
    auto word0 = this->split(delimiter);
    std::vector<text> word;
    for(int i=0;i<word0.size();i++){
        if( word0[i].Contain(equalSign) ){
            if( word0[i].left(1).ContentEqualCaseInsensitive(equalSign) ){
                // bugged
            }else if( word0[i].right(1).ContentEqualCaseInsensitive(equalSign) ){
                if( i + 1 < word0.size() && !word0[i+1].Contain(equalSign) ){
                    word.push_back(word0[i] + word0[i+1]);
                    i += 1;
                }else{
                    // bugged
                }
            }else{
                word.push_back(word0[i]);
            }
        }else if( i + 2 < word0.size() &&
            word0[i+1].ContentEqualCaseInsensitive(equalSign) &&
            !word0[i+2].Contain(equalSign) ){
            word.push_back(word0[i] + word0[i+1] + word0[i+2]);
            i += 2;
        }else if( i + 1 < word0.size() && 
            word0[i+1].Contain(equalSign) ){
            word.push_back(word0[i] + word0[i+1]);
            i += 1;
        }
    }
    for(int i=0;i<word.size();i++){
        auto word1 = word[i].split(equalSign);
        if( word1.size() == 2 ){
            result.emplace_back(word1[0], word1[1]);
        }
    }
    return result;
}
int text::numberOfAppearances(text word) const{
    int count = 0;
    size_t pos = 0;
    while( true ){
        pos = data.find(word.string(), pos);
        if( pos == std::string::npos ) break;
        count++;
        pos += word.size();
    }
    return count;
}

const char* text::c_str() const{
    return data.c_str();
};

const std::string text::string() const{
    return data;
};

text text::FromTextArray(std::vector<text> *word, int first, int last){
    text result;
    if( last > word->size() ) last = word->size();
    for(int i=first;i<last;i++){
        if( i != last - 1 ){
            result += (*word)[i] + " ";
        }else{
            result += (*word)[i];
        }
    }
    return result;
};

text text::FromDoubleArray(std::vector<double> *word, text format, int first, int last){
    text result;
    if( last > word->size() ) last = word->size();
    for(int i=first;i<last;i++){
        result += text(format.c_str(), (*word)[i]) + " ";
    }
    return result;
};

text text::FromIntArray(std::vector<int> *word, text format, int first, int last){
    text result;
    if( last > word->size() ) last = word->size();
    for(int i=first;i<last;i++){
        result += text(format.c_str(), (*word)[i]) + " ";
    }
    return result;
};

text text::FromByte2Readable(double value){
    double value2 = value;
    const text unit[5] = { "B", "K", "M", "G", "T" };
    int iunit = 0;
    while( value2 > 1024.0 ){
        value2 /= 1024.0;
        iunit++;
    }
    char buffer[256];
    sprintf(buffer, "%13.7f %s", value2, unit[iunit].c_str());
    return buffer;
};

text text::FromMilliSecond2Readable(double value){
    double time = value;
    text unit = "ms";

    if( time > 1000.0 ){
        time = time / 1000.0;
        unit = "s";
        if( time > 60.0 ){
            time = time / 60.0;
            unit = "m";
            if( time > 60.0 ){
                time = time / 60.0;
                unit = "h";
                if( time > 24.0 ){
                    time = time / 24.0;
                    unit = "d";
                }
            }
        }
    }

    char buffer[256];
    sprintf(buffer, "%6.2f %3s", time, unit.c_str());
    return buffer;
};

text text::FromSecond2Readable(double value){
    double time = value;
    text unit = "s";

    if( time > 60.0 ){
        time = time / 60.0;
        unit = "m";
        if( time > 60.0 ){
            time = time / 60.0;
            unit = "h";
            if( time > 24.0 ){
                time = time / 24.0;
                unit = "d";
            }
        }
    }

    char buffer[256];
    sprintf(buffer, "%6.2f %3s", time, unit.c_str());
    return buffer;
};

text text::GetLine(std::ifstream *fin){
    std::string buffer;
    getline(*fin, buffer);
    return buffer;
};

text text::formatC(const char *format, ...){
    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
    return buffer;
}

text operator + (const text& lhs, const text& rhs){
    return lhs.string() + rhs.string();
};

bool operator < (const text& lhs, const text& rhs){
    return lhs.string() < rhs.string();
};

bool operator > (const text& lhs, const text& rhs){
    return lhs.string() > rhs.string();
};

bool operator == (const text& lhs, const text& rhs){
    return lhs.string() == rhs.string();
};

bool operator != (const text& lhs, const text& rhs){
    return lhs.string() != rhs.string();
};

text operator * (const text& lhs, const int& rhs){
    text result;
    for(int i=0;i<rhs;i++){
        result += lhs;
    }
    return result;
}

std::ostream& operator<<(std::ostream& lhs, const text& rhs){
    lhs << rhs.string();
    return lhs;
}