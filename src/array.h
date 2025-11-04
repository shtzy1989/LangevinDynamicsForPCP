#ifndef MD_ARRAY_H
#define MD_ARRAY_H

#include "cppheader.h"
#include <vector>
#include <algorithm>

template<typename T>
class Array{
private:
    std::vector<T> m_data;
public:
    // constructors
    Array(){
    };
    Array(size_t size) : m_data(size){ 
    };
    Array(size_t size, const T& value) : m_data(size, value){
    };
    Array<T>(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last) : m_data(first, last){
    };
    Array(const std::vector<T>& x) : m_data(x){
    };
    // copy
    Array<T>& operator= (const Array<T>& x){
        m_data = x.m_data;
        return *this;
    };
    Array<T>& operator= (const std::vector<T>& x){
        m_data = x;
        return *this;
    };
public:
    typename std::vector<T>::iterator begin(){
        return m_data.begin();
    };
    typename std::vector<T>::iterator end(){
        return m_data.end();
    };
    typename std::vector<T>::reverse_iterator rbegin(){
        return m_data.rbegin();
    };
    typename std::vector<T>::reverse_iterator rend(){
        return m_data.rend();
    };
    typename std::vector<T>::iterator cbegin() const{
        return m_data.cbegin();
    };
    typename std::vector<T>::iterator cend() const{
        return m_data.cend();
    };
    typename std::vector<T>::reverse_iterator crbegin() const{
        return m_data.crbegin();
    };
    typename std::vector<T>::reverse_iterator crend() const{
        return m_data.crend();
    };
public:
    size_t size() const{
        return m_data.size();
    };
    size_t max_size() const{
        return m_data.max_size();
    };
    void resize(size_t n){
        m_data.resize(n);
    };
    void resize(size_t n, T value){
        m_data.resize(n, value);
    };
    size_t capacity() const{
        return m_data.capacity();
    };
    bool empty() const{
        return m_data.empty();
    };
    void reserve(size_t n){
        m_data.reserve(n);
    };
    void shrink_to_fit(){
        m_data.shrink_to_fit();
    };
public:
    T& operator[] (size_t n){
        return m_data[n];
    };
    const T& operator[] (size_t n) const{
        return m_data[n];
    };
    T& at(size_t n){
        return m_data.at(n);
    };
    const T& at(size_t n) const{
        return m_data.at(n);
    };
    T& front(){
        return m_data.front();
    };
    const T& front() const{
        return m_data.front();
    };
    T& back(){
        return m_data.back();
    };
    const T& back() const{
        return m_data.back();
    };
    T* data() noexcept{
        return m_data.data();
    };
    const T* data() const noexcept{
        return m_data.data();
    };
public:
    void assign(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        m_data.assign(first, last);
    };
    void assign(size_t n, const T& value){
        m_data.assign(n, value);
    };
    void push_back(const T& value){
        m_data.push_back(value);
    };
    void pop_back(){
        m_data.pop_back();
    };
    typename std::vector<T>::iterator insert(typename std::vector<T>::iterator position, const T& value){
        return m_data.insert(position, value);
    };
    void insert(typename std::vector<T>::iterator position, size_t n, const T& value){
        m_data.insert(position, n, value);
    };
    void insert(typename std::vector<T>::iterator position, typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        m_data.insert(position, first, last);
    };
    typename std::vector<T>::iterator erase(typename std::vector<T>::iterator position){
        return m_data.erase(position);
    };
    typename std::vector<T>::iterator erase(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        return m_data.erase(first, last);
    };
    void swap(std::vector<T>& x){
        m_data.swap(x);
    };
    void swap(Array<T>& x){
        m_data.swap(x.m_data);
    };
    void clear(){
        m_data.clear();
    };
    template <class... Args> typename std::vector<T>::iterator emplace(const typename std::vector<T>::iterator position, Args&&... args){
        return m_data.emplace(position, args...);
    };
    template <class... Args> void emplace_back(Args&&... args){
        m_data.emplace_back(args...);
    };
public:
    std::vector<T> ToVector(){
        return m_data;
    };
    std::vector<T> ToVector() const{
        return m_data;
    };
    std::vector<T>* ToVectorAddress(){
        return &m_data;
    };
    bool Assign(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        try{
            m_data.assign(first, last);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Assign(size_t n, const T& value){
        try{
            m_data.assign(n, value);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Resize(int n){
        try{
            m_data.resize(n);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Resize(int n, const T& value){
        try{
            m_data.resize(n, value);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Reserve(int n){
        try{
            m_data.reserve(n);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Add(const T& value){
        try{
            m_data.push_back(value);
        }catch(...){
            return false;
        }
        return true;
    };
    template <class... Args> bool AddArgument(Args&&... args){
        try{
            m_data.emplace_back(args...);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Insert(typename std::vector<T>::iterator position, const T& value){
        try{
            m_data.insert(position, value);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Insert(typename std::vector<T>::iterator position, size_t n, const T& value){
        try{
            m_data.insert(position, n, value);
        }catch(...){
            return false;
        }
        return true;
    };
    bool Insert(typename std::vector<T>::iterator position, typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        try{
            m_data.insert(position, first, last);
        }catch(...){
            return false;
        }
        return true;
    };
    
    T Sum(){
        return Sum(m_data.begin(), m_data.end());
    };
    T Sum(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        T value = 0;
        for(auto i=first;i!=last;i++){
            value = value + *i;
        }
        return value;
    };
    T Average(){
        return Average(m_data.begin(), m_data.end());
    };
    T Average(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        T value = 0;
        int count = 0;
        for(auto i=first;i!=last;i++){
            value = value + *i;
            count++;
        }
        value /= count;
        return value;
    }; 
    T RMSD(){
        return RMSD(m_data.begin(), m_data.end());
    };
    T RMSD(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        T avg = 0;
        T avg2 = 0;
        for(auto i=first;i!=last;i++){
            avg = avg + *i;
            avg2 = avg2 + *i * *i;
        }
        avg /= m_data.size();
        avg2 /= m_data.size();

        return sqrt(avg2 - avg * avg);
    };
    T Max(){
        return Max(m_data.begin(), m_data.end());
    };
    T Max(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        T value = *first;
        for(auto i=first;i!=last;i++){
            value = *i > value ? *i : value;
        }
        return value;
    };
    T Min(){
        return Min(m_data.begin(), m_data.end());
    };
    T Min(typename std::vector<T>::iterator first, typename std::vector<T>::iterator last){
        T value = *first;
        for(auto i=first;i!=last;i++){
            value = *i < value ? *i : value;
        }
        return value;
    };
    void Print(FILE* file = stderr, int column = 1, const char *format = "%13.7f "){
        Print(file, m_data.begin(), m_data.end(), column, format);
    };
    void Print(FILE* file, typename std::vector<T>::iterator first, typename std::vector<T>::iterator last, int column = 1, const char *format = "%13.7f "){
        if( column <= 0 ) column = 1;
        int count = 0;
        bool returned = false;
        for(auto i=first;i!=last;i++){
            fprintf(file, format, *i);
            returned = false;
            count++;
            if( count % column == 0 ){
                fprintf(file, "\n");
                returned = true;
            }
        }
        if( !returned ){
            fprintf(file, "\n");
        }
    };
    bool InitializeArray(const T& first, const T& last, const T& interval){
        if( interval <= 0 ) return false;
        int size = ceil(double(last - first) / interval);
        bool hr = Resize(size);
        if( !hr ){
            return false;
        }else{
            int position = 0;
            for(auto i=first;i<last;i+=interval){
                m_data[position++] = i;
            }
            return true;
        }
    };
    static Array<T> InitializeArrayS(const T& first, const T& last, const T& interval){
        Array<T> result;
        if( interval <= 0 ) return result;
        int size = ceil(double(last - first) / interval);
        bool hr = result.Resize(size);
        if( !hr ){
            return result;
        }else{
            int position = 0;
            for(auto i=first;i<last;i+=interval){
                result.m_data[position++] = i;
            }
            return result;
        }
    };
    Array<T> SubArray(size_t first, size_t last){
        Array<T> result;
        result.Assign(m_data.begin() + first, m_data.begin() + last);
        return result;
    };
    void Sort(){
        std::sort(m_data.begin(), m_data.end());
    };
    bool Contain(const T& value){
        bool in = false; 
        for(int i=0;i<m_data.size();i++){
            if( m_data[i] == value ){
                in = true; 
                break;
            }
        }
        return in;
    };
    int FindIndex(const T& value){
        int in = -1;
        for(int i=0;i<m_data.size();i++){
            if( m_data[i] == value ){
                in = i;
                break;
            }
        }
        return in;
    };
};

template<typename T> Array<T> operator + (Array<T>& lhs, Array<T>& rhs){
    Array<T> result;
    result.Insert(result.end(), lhs.begin(), lhs.end());
    result.Insert(result.end(), rhs.begin(), rhs.end());
    return result;
};

template<typename T> bool operator ==(const Array<T>& a, const Array<T>& b){
    bool equal = a.size() == b.size();
    if( !equal ){
        return false;
    }else{
        for(int i=0;i<a.size();i++){
            if( a[i] != b[i] ){
                equal = false;
                break;
            }
        }
        return equal;
    }
};

template<typename T> bool operator !=(const Array<T>& a, const Array<T>& b){
    return !(a == b);
};

template<typename T> int BinarySearch(Array<T> *data, T value) {
    if( data->size() == 0 ){
        return -1;
    }    
    int left = 0;
    int right = data->size() - 1;
    int pos;

    while (left <= right) {
        int mid = (left + right) / 2;
        T midValue = (*data)[mid];
        if (midValue < value) {
            left = mid + 1;
            pos = mid + 1;
        }
        else if (midValue > value) {
            right = mid - 1;
            pos = mid;
        }
        else {
            return mid;
        }
    }
    return -1;
};

template<typename T> int BinarySearchReturnPositionBefore(Array<T> *data, T value) {
    if( data->size() == 0 ){
        return 0;
    }    
    int left = 0;
    int right = data->size() - 1;
    int pos;

    while (left <= right) {
        int mid = (left + right) / 2;
        T midValue = (*data)[mid];
        if (midValue < value) {
            left = mid + 1;
            pos = mid + 1;
        }
        else if (midValue > value) {
            right = mid - 1;
            pos = mid;
        }
        else {
            return mid;
        }
    }
    return pos;
};

#endif
