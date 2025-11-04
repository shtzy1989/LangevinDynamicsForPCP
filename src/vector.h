#ifndef MD_MATH_VECTOR_H
#define MD_MATH_VECTOR_H

#include "mathdef.h"
#include "cppheader.h"

template<typename T>
union GeneralVector2D{
	T m[2];
	struct{
		T x, y;
	};
	GeneralVector2D(){};
	GeneralVector2D(T _x, T _y);
	GeneralVector2D(T _x);

	GeneralVector2D<T>& operator =(const GeneralVector2D& rhs);
	GeneralVector2D<T>& operator =(const T *rhs);
	T& operator [] (int index);

	T length() const;
	T length2() const;

	GeneralVector2D<T> normalize();

};

template<typename T> GeneralVector2D<T>::GeneralVector2D(T _x, T _y){
	x = _x;
	y = _y;
};

template<typename T> GeneralVector2D<T>::GeneralVector2D(T _x){
	x = _x;
	y = _x;
};

template<typename T> GeneralVector2D<T>& GeneralVector2D<T>::operator =(const GeneralVector2D& rhs){
	x = rhs.x;
	y = rhs.y;
	return *this;
};

template<typename T> GeneralVector2D<T>& GeneralVector2D<T>::operator =(const T *rhs){
	x = rhs[0];
	y = rhs[1];
	return *this;
};

template<typename T> T& GeneralVector2D<T>::operator [] (int index){
	return m[index];
};

template<typename T> T GeneralVector2D<T>::length() const{
	return std::sqrt(x * x + y * y);
};

template<typename T> T GeneralVector2D<T>::length2() const{
	return x * x + y * y;
};

template<typename T> GeneralVector2D<T> GeneralVector2D<T>::normalize(){
	GeneralVector2D<T> result;
	double dnorm = length();
	result.x = x / dnorm;
	result.y = y / dnorm;
	return result;
};

template<typename T> bool operator == (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	return lhs.x == rhs.x && lhs.y == rhs.y;
};

template<typename T> bool operator != (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	return !(lhs == rhs);
};

template<typename T> GeneralVector2D<T> operator + (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	GeneralVector2D<T> result;
	result.x = lhs.x + rhs.x;
	result.y = lhs.y + rhs.y;
	return result;
};

template<typename T> GeneralVector2D<T> operator - (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	GeneralVector2D<T> result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	return result;
};

template<typename T> T operator * (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	return lhs.x * rhs.x + lhs.y * rhs.y;
};

template<typename T> GeneralVector2D<T> operator * (const GeneralVector2D<T>& lhs, const T& rhs){
	GeneralVector2D<T> result;
	result.x = lhs.x * rhs;
	result.y = lhs.y * rhs;
	return result;
};

template<typename T> GeneralVector2D<T> VecCrossProduct (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	GeneralVector2D<T> result;
	result.x = lhs.y * rhs.z - lhs.z * rhs.y;
	result.y = lhs.z * rhs.x - lhs.x * rhs.z;
	return result;
};

template<typename T> GeneralVector2D<T> operator / (const GeneralVector2D<T>& lhs, const T& rhs){
	GeneralVector2D<T> result;
	result.x = lhs.x / rhs;
	result.y = lhs.y / rhs;
	return result;
};

template<typename T> T VecCrossDegree (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	T degree;
	T n1 = lhs.length2();
	T n2 = rhs.length2();
	if( n1 == 0 || n2 == 0 ) return 0.0;
	T v = lhs * rhs / std::sqrt(n1 * n2);
	if( v > 1.0 ) v = 1.0;
	if( v < -1.0 ) v = -1.0;
	degree = std::acos(v);
	return degree;
};

template<typename T> T VecDistance (const GeneralVector2D<T>& lhs, const GeneralVector2D<T>& rhs){
	GeneralVector2D<T> sub = lhs - rhs;
	return sub.Norm();
};

template<typename T> T VecDihedral (const GeneralVector2D<T>& v1, const GeneralVector2D<T>& v2, const GeneralVector2D<T>& v3, const GeneralVector2D<T>& v4){
    GeneralVector2D<T> l1, l2, l3, n;
    l1 = v1 - v2;
    l2 = v2 - v3;
    l3 = v3 - v4;
    l1 = VecCrossProduct(l1, l2);
    l3 = VecCrossProduct(l2, l3);
    n = VecCrossProduct(l1, l3);
    if( n * l2 < 0 )
        return VecCrossDegree(l1, l3);
    else
        return -VecCrossDegree(l1, l3);
};

typedef GeneralVector2D<double> DVec2D;
typedef GeneralVector2D<int> IVec2D;
typedef GeneralVector2D<float> FVec2D;
typedef GeneralVector2D<Real> RVec2D;

bool operator < (const IVec2D& lhs, const IVec2D& rhs);

// GeneralVector ========================================================================

template<typename T>
union GeneralVector{
	T m[3];
	struct{
		T x, y, z;
	};
	GeneralVector(){};
	GeneralVector(T _x, T _y, T _z);
	GeneralVector(T _x);

	GeneralVector<T>& operator =(const GeneralVector& rhs);
	GeneralVector<T>& operator =(const T *rhs);
	T& operator [] (int index);

	T length() const;
	T length2() const;

	GeneralVector<T> normalize();

};

template<typename T> GeneralVector<T>::GeneralVector(T _x, T _y, T _z){
	x = _x;
	y = _y;
	z = _z;
};

template<typename T> GeneralVector<T>::GeneralVector(T _x){
	x = _x;
	y = _x;
	z = _x;
};

template<typename T> GeneralVector<T>& GeneralVector<T>::operator =(const GeneralVector& rhs){
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	return *this;
};

template<typename T> GeneralVector<T>& GeneralVector<T>::operator =(const T *rhs){
	x = rhs[0];
	y = rhs[1];
	z = rhs[2];
	return *this;
};

template<typename T> T& GeneralVector<T>::operator [] (int index){
	return m[index];
};

template<typename T> T GeneralVector<T>::length() const{
	return std::sqrt(x * x + y * y + z * z);
};

template<typename T> T GeneralVector<T>::length2() const{
	return x * x + y * y + z * z;
};

template<typename T> GeneralVector<T> GeneralVector<T>::normalize(){
	GeneralVector<T> result;
	double dnorm = length();
	result.x = x / dnorm;
	result.y = y / dnorm;
	result.z = z / dnorm;
	return result;
};

template<typename T> bool operator == (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
};

template<typename T> bool operator != (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	return !(lhs == rhs);
};

template<typename T> GeneralVector<T> operator + (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	GeneralVector<T> result;
	result.x = lhs.x + rhs.x;
	result.y = lhs.y + rhs.y;
	result.z = lhs.z + rhs.z;
	return result;
};

template<typename T> GeneralVector<T> operator - (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	GeneralVector<T> result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	result.z = lhs.z - rhs.z;
	return result;
};

template<typename T> T operator * (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
};

template<typename T> GeneralVector<T> operator * (const GeneralVector<T>& lhs, const T& rhs){
	GeneralVector<T> result;
	result.x = lhs.x * rhs;
	result.y = lhs.y * rhs;
	result.z = lhs.z * rhs;
	return result;
};

template<typename T> GeneralVector<T> VecCrossProduct (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	GeneralVector<T> result;
	result.x = lhs.y * rhs.z - lhs.z * rhs.y;
	result.y = lhs.z * rhs.x - lhs.x * rhs.z;
	result.z = lhs.x * rhs.y - lhs.y * rhs.x;
	return result;
};

template<typename T> GeneralVector<T> operator / (const GeneralVector<T>& lhs, const T& rhs){
	GeneralVector<T> result;
	result.x = lhs.x / rhs;
	result.y = lhs.y / rhs;
	result.z = lhs.z / rhs;
	return result;
};

template<typename T> T VecCrossDegree (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	T degree;
	T n1 = lhs.length2();
	T n2 = rhs.length2();
	if( n1 == 0 || n2 == 0 ) return 0.0;
	T v = lhs * rhs / std::sqrt(n1 * n2);
	if( v > 1.0 ) v = 1.0;
	if( v < -1.0 ) v = -1.0;
	degree = std::acos(v);
	return degree;
};

template<typename T> T VecDistance (const GeneralVector<T>& lhs, const GeneralVector<T>& rhs){
	GeneralVector<T> sub = lhs - rhs;
	return sub.Norm();
};

template<typename T> T VecDihedral (const GeneralVector<T>& v1, const GeneralVector<T>& v2, const GeneralVector<T>& v3, const GeneralVector<T>& v4){
    GeneralVector<T> l1, l2, l3, n;
    l1 = v1 - v2;
    l2 = v2 - v3;
    l3 = v3 - v4;
    l1 = VecCrossProduct(l1, l2);
    l3 = VecCrossProduct(l2, l3);
    n = VecCrossProduct(l1, l3);
    if( n * l2 < 0 )
        return VecCrossDegree(l1, l3);
    else
        return -VecCrossDegree(l1, l3);
};

typedef GeneralVector<double> DVec;
typedef GeneralVector<int> IVec;
typedef GeneralVector<float> FVec;
typedef GeneralVector<Real> RVec;

RVec operator * (const RVec& lhs, const IVec& rhs);
RVec operator * (const IVec& lhs, const RVec& rhs);

bool operator < (const IVec& lhs, const IVec& rhs);

// ===========================
// 9 dimension vector
// ===========================

template<typename T>
union GeneralVector9{
	T m[9];
	struct{
		T xx, xy, xz, yx, yy, yz, zx, zy, zz;
	};
	GeneralVector9(){};
	GeneralVector9(T _xx, T _xy, T _xz, T _yx, T _yy, T _yz, T _zx, T _zy, T _zz);

	GeneralVector9<T>& operator =(const GeneralVector9& rhs);
	GeneralVector9<T>& operator =(const T *rhs);
	T& operator [] (int index);
	T& Diagonal(int index);
	GeneralVector<T> DiagonalToVec();

	void Print(const char *format = "%13.7f ", FILE *fout = stderr);
};

template<typename T> GeneralVector9<T>::GeneralVector9(T _xx, T _xy, T _xz, T _yx, T _yy, T _yz, T _zx, T _zy, T _zz){
	xx = _xx;
	xy = _xy;
	xz = _xz;
	yx = _yx;
	yy = _yy;
	yz = _yz;
	zx = _zx;
	zy = _zy;
	zz = _zz;
};

template<typename T> GeneralVector9<T>& GeneralVector9<T>::operator =(const GeneralVector9& rhs){
	xx = rhs.xx;
	xy = rhs.xy;
	xz = rhs.xz;
	yx = rhs.yx;
	yy = rhs.yy;
	yz = rhs.yz;
	zx = rhs.zx;
	zy = rhs.zy;
	zz = rhs.zz;
	return *this;
};

template<typename T> GeneralVector9<T>& GeneralVector9<T>::operator =(const T *rhs){
	xx = rhs[0];
	xy = rhs[1];
	xz = rhs[2];
	yx = rhs[3];
	yy = rhs[4];
	yz = rhs[5];
	zx = rhs[6];
	zy = rhs[7];
	zz = rhs[8];
	return *this;
};

template<typename T> T& GeneralVector9<T>::operator [] (int index){
	return m[index];
};

template<typename T> T& GeneralVector9<T>::Diagonal(int index){
	return m[index*3+index];
};

template<typename T> 
GeneralVector<T> GeneralVector9<T>::DiagonalToVec(){
	return RVec(xx, yy, zz);
};
template<typename T> void GeneralVector9<T>::Print(const char *format, FILE *fout){
	fprintf(fout, format, xx);
	fprintf(fout, format, xy);
	fprintf(fout, format, xz);
	fprintf(fout, "\n");
	fprintf(fout, format, yx);
	fprintf(fout, format, yy);
	fprintf(fout, format, yz);
	fprintf(fout, "\n");
	fprintf(fout, format, zx);
	fprintf(fout, format, zy);
	fprintf(fout, format, zz);
	fprintf(fout, "\n");
}

typedef GeneralVector9<double> DVec9;
typedef GeneralVector9<int> IVec9;
typedef GeneralVector9<float> FVec9;
typedef GeneralVector9<Real> RVec9;


template<typename T, size_t dimension>
class GeneralDimensionVector{
private:
    std::vector<T> m_Data;
public:
    GeneralDimensionVector(){
        m_Data.resize(dimension);
		// memset(m_Data, 0, sizeof(T) * dimension);
    };
    GeneralDimensionVector(T* value){
        m_Data.resize(dimension);
        for(int i=0;i<dimension;i++){
            m_Data[i] = value[i];
        }
    }
    GeneralDimensionVector(std::vector<T> value){
        m_Data.resize(dimension);
        for(int i=0;i<dimension;i++){
            m_Data[i] = value[i];
        }
    }
    ~GeneralDimensionVector(){
        //SAFE_DELETE_ARRAY(m_Data);
    }
public:
    size_t Dimension() const{
        return dimension;
    }
public:
    GeneralDimensionVector& operator = (const GeneralDimensionVector& rhs){
        for(int i=0;i<dimension;i++){
            m_Data[i] = rhs.m_Data[i];
        }
        return *this;
    };
    GeneralDimensionVector& operator += (const GeneralDimensionVector& rhs){
        for(int i=0;i<dimension;i++){
            m_Data[i] += rhs.m_Data[i];
        }
        return *this;
    }
    GeneralDimensionVector& operator -= (const GeneralDimensionVector& rhs){
        for(int i=0;i<dimension;i++){
            m_Data[i] -= rhs.m_Data[i];
        }
        return *this;
    }
    const T& operator [](size_t p) const{
        return m_Data[p];
    }
    T& operator [](size_t p){
        return m_Data[p];
    }
    bool operator < (const GeneralDimensionVector rhs) const{
        for(int i=0;i<dimension;i++){
            if( m_Data[i] < rhs.m_Data[i] ) return true;
            else if( m_Data[i] > rhs.m_Data[i] ) return false;
        }
        return false;
    }
    bool operator > (const GeneralDimensionVector rhs) const{
        for(int i=0;i<dimension;i++){
            if( m_Data[i] > rhs.m_Data[i] ) return true;
            else if( m_Data[i] < rhs.m_Data[i] ) return false;
        }
        return false;
    }
    T length2(){
        T value = 0.0;
        for(int i=0;i<dimension;i++){
            value += m_Data[i] * m_Data[i];
        }
        return value;
    }
    T length(){
        return sqrt(length2());
    }
    GeneralDimensionVector normalizeSelf(){
        T length = length();
        if( length == 0 ){
            return *this;
        }else{
            for(int i=0;i<dimension;i++){
                m_Data[i] /= length;
            }
            return *this;
        }
    }
    GeneralDimensionVector normalize(){
        GeneralDimensionVector result = *this;
        T length = result.length();
        if( length == 0 ){
            return result;
        }else{
            for(int i=0;i<dimension;i++){
                result.m_Data[i] /= length;
            }
            return result;
        }
    }
    GeneralDimensionVector scaleSelf(T value){
        for(int i=0;i<dimension;i++){
            m_Data[i] *= value;
        }
        return *this;
    }
    GeneralDimensionVector scale(T value){
        GeneralDimensionVector result = *this;
        for(int i=0;i<dimension;i++){
            result.m_Data[i] *= value;
        }
        return result;
    }
    T sum(){
        T value = 0.0;
        for(int i=0;i<dimension;i++){
            value += m_Data[i];
        }
        return value;
    }
    GeneralDimensionVector zeroSelf(){
        for(int i=0;i<dimension;i++){
            m_Data[i] = 0.0;
        }
        return *this;
    }
	void Print(FILE* fout = stderr, const char* channelFormat = "%13.7f "){
		for(int i=0;i<dimension;i++){
            fprintf(fout, channelFormat, m_Data[i]);
        }
		fprintf(fout, "\n");
	}
	void PrintRich(FILE* fout = stderr, const char* overFormat = "%s\n", const char* channelFormat = "%13.7f "){
		char buffer[1024] = "";
		for(int i=0;i<dimension;i++){
			char buffer2[256];
            sprintf(buffer2, channelFormat, m_Data[i]);
			strcat(buffer, buffer2);
        }
		fprintf(fout, overFormat, buffer);
	}
};

template<typename T, size_t dimension>
bool operator ==(GeneralDimensionVector<T, dimension> a, GeneralDimensionVector<T, dimension> b){
	bool bResult = true;
	for(int i=0;i<dimension;i++){
		bResult = bResult && a[i] == b[i];
	}
	return bResult;
};

template<typename T, size_t dimension>
bool operator !=(GeneralDimensionVector<T, dimension> a, GeneralDimensionVector<T, dimension> b){
	return !(a == b);
};

template<typename T, size_t dimension>
GeneralDimensionVector<T, dimension> operator + (const GeneralDimensionVector<T, dimension>& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	GeneralDimensionVector<T, dimension> result;
	for(int i=0;i<dimension;i++){
		result[i] = lhs[i] + rhs[i];
	}
	return result;
};

template<typename T, size_t dimension>
GeneralDimensionVector<T, dimension> operator - (const GeneralDimensionVector<T, dimension>& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	GeneralDimensionVector<T, dimension> result;
	for(int i=0;i<dimension;i++){
		result[i] = lhs[i] - rhs[i];
	}
	return result;
};

template<typename T, size_t dimension>
T operator * (const GeneralDimensionVector<T, dimension>& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	T value = 0;
	for(int i=0;i<dimension;i++){
		value += lhs[i] * rhs[i];
	}
	return value;
};

template<typename T, size_t dimension>
GeneralDimensionVector<T, dimension> operator * (const GeneralDimensionVector<T, dimension>& lhs, const T& rhs){
	GeneralDimensionVector<T, dimension> result;
	for(int i=0;i<dimension;i++){
		result[i] = lhs[i] * rhs;
	}
	return result;
};

template<typename T, size_t dimension>
GeneralDimensionVector<T, dimension> operator * (const T& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	GeneralDimensionVector<T, dimension> result;
	for(int i=0;i<dimension;i++){
		result[i] = lhs * rhs[i];
	}
	return result;
}

template<typename T, size_t dimension>
GeneralDimensionVector<T, dimension> operator / (const GeneralDimensionVector<T, dimension>& lhs, const T& rhs){
	GeneralDimensionVector<T, dimension> result;
	for(int i=0;i<dimension;i++){
		result[i] = lhs[i] / rhs;
	}
	return result;
}

template<typename T, size_t dimension>
T GeneralDimensionVectorDistance(const GeneralDimensionVector<T, dimension>& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	GeneralDimensionVector<T, dimension> subtract = lhs - rhs;
	return subtract.length();
}

template<typename T, size_t dimension>
T GeneralDimensionVectorDistance2(const GeneralDimensionVector<T, dimension>& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	GeneralDimensionVector<T, dimension> subtract = lhs - rhs;
	return subtract.length2();
}

template<typename T, size_t dimension>
GeneralDimensionVector<T, dimension> GeneralDimensionVectorChannelMultiply (const GeneralDimensionVector<T, dimension>& lhs, const GeneralDimensionVector<T, dimension>& rhs){
	GeneralDimensionVector<T, dimension> result;
	for(int i=0;i<dimension;i++){
		result[i] = lhs[i] * rhs[i];
	}
	return result;
}

typedef GeneralDimensionVector<double, 2> GDDoubleVec2;
typedef GeneralDimensionVector<double, 3> GDDoubleVec3;

#endif
