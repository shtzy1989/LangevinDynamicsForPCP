#include "vector.h"


bool operator < (const IVec2D& lhs, const IVec2D& rhs){
	return (lhs.x == rhs.x && lhs.y < rhs.y) || lhs.x < rhs.x;
};

// GeneralVector ==================================================================
RVec operator * (const RVec& lhs, const IVec& rhs) {
	return RVec(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
};

RVec operator * (const IVec& lhs, const RVec& rhs) {
	return RVec(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
};

bool operator < (const IVec& lhs, const IVec& rhs){
	return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z < rhs.z) || (lhs.x == rhs.x && lhs.y < rhs.y) || lhs.x < rhs.x;
};