#include "simd.h"
#include "cppheader.h"

void printmm512d(char *tag, __m512d* value){
    fprintf(stderr, "%20s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
        tag,
        ((double*)value)[0], ((double*)value)[1], ((double*)value)[2], ((double*)value)[3],
        ((double*)value)[4], ((double*)value)[5], ((double*)value)[6], ((double*)value)[7]);
};

void print512k(int k){
    fprintf(stderr, "%20s %8d %8d %8d %8d %8d %8d %8d %8d\n",
        "index", k, k+1, k+2, k+3, k+4, k+5, k+6, k+7);
};


void printmm256d(char *tag, __m256d* value){
    fprintf(stderr, "%20s %13.7f %13.7f %13.7f %13.7f\n", tag, ((double*)value)[0], ((double*)value)[1], ((double*)value)[2], ((double*)value)[3]);
};

void print256k(int k){
    fprintf(stderr, "%20s %13d %13d %13d %13d\n", "index", k, k+1, k+2, k+3);
};

#ifdef SIMDLENGTH8
__m512d _mmSIMD_fmadd_pd(__m512d __A, __m512d __B, __m512d __C){
    return _mm512_fmadd_pd(__A, __B, __C);
}
#else
// this does not work because _mm256_cmp_pd requires an immediate _P
// __mmask8 _mm256_cmp_pd_mask_avx256(__m256d _X, __m256d _Y, const int _P){
//     __m256d _mask256 = _mm256_cmp_pd(_X, _Y, _P);
//     __mmask8 _mask = _mm256_movemask_pd(_mask256);
//     return _mask;
// };
__m256d _mm256_mask_blend_pd_avx256(__mmask8 _U, __m256d _A, __m256d _W){
    __m256d _mask = _mm256_set_pd(
        (_U & 8) ? -1.0 : 0.0,
        (_U & 4) ? -1.0 : 0.0,
        (_U & 2) ? -1.0 : 0.0,
        (_U & 1) ? -1.0 : 0.0
    );
    return _mm256_blendv_pd(_A, _W, _mask);
};
__mmask8 _mmask8_and_SIMD(__mmask8 _a, __mmask8 _b){
    return _a & _b;
};
__mmask8 _mmask8_andn_SIMD(__mmask8 _a, __mmask8 _b){
    return ~_a & _b;
}
#endif


