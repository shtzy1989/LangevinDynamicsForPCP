#ifndef MD_SIMD_H
#define MD_SIMD_H

#include "cppheader.h"

#include "simd512.h"

#include <immintrin.h>

#    define GMX_DOUBLE_NEGZERO  (-0.0)
#define GMX_SIMD_RCP_BITS                      11
#define GMX_SIMD_ACCURACY_BITS_DOUBLE 44
// simd/impl_x86_avx_256/impl_x86_avx_256_simd_double.h
// simd/simd_math.h

inline static void printBinary8(const char *info, char a){
    int i;
    fprintf(stderr, "%20s 0b", info);
    for (int i = 7; i >= 0; i--) {
        fprintf(stderr, "%d", (a >> i) & 1);
    }
    fprintf(stderr, "\n");
};
inline static void _mmask8_store_256(bool* output, __mmask8 _input){
    output[0] = _input & 0b00000001;
    output[1] = _input & 0b00000010;
    output[2] = _input & 0b00000100;
    output[3] = _input & 0b00001000;
};

inline static void _mmask8_store_512(bool* output, __mmask8 _input){
    output[0] = _input & 0b00000001;
    output[1] = _input & 0b00000010;
    output[2] = _input & 0b00000100;
    output[3] = _input & 0b00001000;
    output[4] = _input & 0b00010000;
    output[5] = _input & 0b00100000;
    output[6] = _input & 0b01000000;
    output[7] = _input & 0b10000000;
};

inline static __mmask8 _mmask8_load_256(bool* input){
    __mmask8 result = 0;
    // if( input[0] ) result = result | 0b00000001;
    // if( input[1] ) result = result | 0b00000010;
    // if( input[2] ) result = result | 0b00000100;
    // if( input[3] ) result = result | 0b00001000;
    result |= input[3]<<3;
    result |= input[2]<<2;
    result |= input[1]<<1;
    result |= input[0];
    return result;
};

inline static __mmask8 _mmask8_load_512(bool* input){
    __mmask8 result = 0;
    // result = result | 0b00000001;
    // result = result | 0b00000010;
    // result = result | 0b00000100;
    // result = result | 0b00001000;
    // result = result | 0b00010000;
    // result = result | 0b00100000;
    // result = result | 0b01000000;
    // result = result | 0b10000000;
    result |= input[7]<<7;
    result |= input[6]<<6;
    result |= input[5]<<5;
    result |= input[4]<<4;
    result |= input[3]<<3;
    result |= input[2]<<2;
    result |= input[1]<<1;
    result |= input[0];
    return result;
};

void printmm512d(char *tag, __m512d* value);

void print512k(int k);


void printmm256d(char *tag, __m256d* value);
void printk256(int k);

static inline __m256d 
fnma(__m256d a, __m256d b, __m256d c)
{
    return _mm256_sub_pd(c, _mm256_mul_pd(a, b));
}

static inline __m256d
fma(__m256d a, __m256d b, __m256d c)
{
    return _mm256_add_pd(_mm256_mul_pd(a, b), c);
}

static inline __m256d
abs(__m256d x)
{
    return _mm256_andnot_pd( _mm256_set1_pd(GMX_DOUBLE_NEGZERO), x );
}

#ifdef SIMDLENGTH8
    const int SIMDWIDTH = 8;
    const int SIMDBYTESIZE = 64;
    typedef __m512d __mSIMDd;
    const auto _mmSIMD_abs_pd = _mm512_abs_pd;
    const auto _mmSIMD_set1_pd = _mm512_set1_pd;
    const auto _mmSIMD_store_pd = _mm512_store_pd;
    const auto _mmSIMD_load_pd = _mm512_load_pd;
    const auto _mmSIMD_add_pd = _mm512_add_pd;
    const auto _mmSIMD_sub_pd = _mm512_sub_pd;
    const auto _mmSIMD_mul_pd = _mm512_mul_pd;
    const auto _mmSIMD_div_pd = _mm512_div_pd;
    const auto _mmSIMD_sqrt_pd = _mm512_sqrt_pd;
#if defined(__INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER) || defined (__INTEL_CLANG_COMPILER)
    // const auto _mmSIMD_log_pd = _mm512_log_pd;
    #define _mmSIMD_log_pd(a) _mm512_log_pd(a)
#endif
    #define _mmSIMD_roundscale_pd(a, b) _mm512_roundscale_pd(a, b)
    #define _mmSIMD_round_pd(a, b) _mm512_roundscale_pd(a, b)
    // const auto _mmSIMD_roundscale_pd = _mm512_roundscale_pd;
    // const auto _mmSIMD_round_pd = _mm512_roundscale_pd;
    const auto _mmSIMD_and_pd = _mm512_and_pd;
    // const auto _mmSIMD_cmp_pd_mask = _mm512_cmp_pd_mask;
    #define _mmSIMD_cmp_pd_mask(a, b, c) _mm512_cmp_pd_mask(a, b, c)
    const auto _mmSIMD_mask_blend_pd = _mm512_mask_blend_pd;
    const auto _mmask8_store_SIMD = _mmask8_store_512;
    const auto _mmask8_load_SIMD = _mmask8_load_512;
    __m512d _mmSIMD_fmadd_pd(__m512d __A, __m512d __B, __m512d __C);
    const auto printmmSIMDd = printmm512d;
    const auto _mmask8_and_SIMD = _kand_mask8;
    const auto _mmask8_andn_SIMD = _kandn_mask8;
#if defined (__INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER) || defined (__INTEL_CLANG_COMPILER)
    #define _mmSIMD_exp_pd(a) _mm512_exp_pd(a)
#endif
    #define _mmSIMD_cos_pd(a) _mm512_cos_pd(a)
    #define _mmSIMD_sin_pd(a) _mm512_sin_pd(a)
#else
    const int SIMDWIDTH = 4;
    const int SIMDBYTESIZE = 32;
    typedef __m256d __mSIMDd;
    #define _mm256_abs_pd(x) _mm256_and_pd((x), _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF)))
    #define _mmSIMD_abs_pd(x) _mm256_abs_pd(x)
    const auto _mmSIMD_set1_pd = _mm256_set1_pd;
    const auto _mmSIMD_store_pd = _mm256_store_pd;
    const auto _mmSIMD_load_pd = _mm256_load_pd;
    const auto _mmSIMD_add_pd = _mm256_add_pd;
    const auto _mmSIMD_sub_pd = _mm256_sub_pd;
    const auto _mmSIMD_mul_pd = _mm256_mul_pd;
    const auto _mmSIMD_div_pd = _mm256_div_pd;
    const auto _mmSIMD_sqrt_pd = _mm256_sqrt_pd;
#if defined(__INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER) || defined (__INTEL_CLANG_COMPILER)
    // const auto _mmSIMD_log_pd = _mm256_log_pd;
    #define _mmSIMD_log_pd(a) _mm256_log_pd(a)
#endif
#ifdef SIMDLENGTH8
    const auto _mmSIMD_roundscale_pd = _mm256_roundscale_pd; // this only exists in avx512
#endif
    #define _mmSIMD_round_pd(a, b) _mm256_round_pd(a, b);
    const auto _mmSIMD_and_pd = _mm256_and_pd;
    // this two exist only in avx512, so defined them
    // const auto _mmSIMD_cmp_pd_mask = _mm256_cmp_pd_mask;
    // const auto _mmSIMD_mask_blend_pd = _mm256_mask_blend_pd; 
    #define _mm256_cmp_pd_mask_avx256(a, b, op) \
        ((unsigned char)_mm256_movemask_pd(_mm256_cmp_pd(a, b, op)))

    __m256d _mm256_mask_blend_pd_avx256(__mmask8 _U, __m256d _A, __m256d _W);
    #define _mmSIMD_cmp_pd_mask(a, b, op) _mm256_cmp_pd_mask_avx256(a, b, op)
    const auto _mmSIMD_mask_blend_pd = _mm256_mask_blend_pd_avx256; 
    
    const auto _mmask8_store_SIMD = _mmask8_store_256;
    const auto _mmask8_load_SIMD = _mmask8_load_256;
    const auto _mmSIMD_fmadd_pd = _mm256_fmadd_pd;
    const auto printmmSIMDd = printmm256d;
    __mmask8 _mmask8_and_SIMD(__mmask8 _a, __mmask8 _b);
    __mmask8 _mmask8_andn_SIMD(__mmask8 _a, __mmask8 _b);
#if defined (__INTEL_COMPILER) || defined (__INTEL_LLVM_COMPILER) || defined (__INTEL_CLANG_COMPILER)
    #define _mmSIMD_exp_pd(a) _mm256_exp_pd(a)
#endif
    #define _mmSIMD_cos_pd(a) _mm256_cos_pd(a)
    #define _mmSIMD_sin_pd(a) _mm256_sin_pd(a)
#endif

#endif

