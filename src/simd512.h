#ifndef SIMD512
#define SIMD512

#include <immintrin.h>

#ifdef SIMDLENGTH8

static inline __mmask16 avx512Int2Mask(short i)
{
    return static_cast<__mmask16>(i);
}

static inline __m512d fma(__m512d a, __m512d b, __m512d c)
{
    return _mm512_fmadd_pd(a, b, c);
}

static inline __m512d ldexp(__m512d value, __m256i exponent)
{
    const __m256i exponentBias = _mm256_set1_epi32(1023);
    __m256i       iExponent    = _mm256_add_epi32(exponent, exponentBias);
    __m512i       iExponent512;

    // Make sure biased argument is not negative
    iExponent = _mm256_max_epi32(iExponent, _mm256_setzero_si256());

    iExponent512 = _mm512_permutexvar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), _mm512_castsi256_si512(iExponent));
    iExponent512 = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), avx512Int2Mask(0xAAAA), iExponent512, 20);
    return _mm512_mul_pd(_mm512_castsi512_pd(iExponent512), value);
}

static inline __m512d
exp(__m512d x)
{
    const __m512d  argscale = _mm512_set1_pd(1.44269504088896340735992468100);
    const __m512d  invargscale0 = _mm512_set1_pd(-0.69314718055966295651160180568695068359375);
    const __m512d  invargscale1 = _mm512_set1_pd(-2.8235290563031577122588448175013436025525412068e-13);
    const __m512d  CE12 = _mm512_set1_pd(2.078375306791423699350304e-09);
    const __m512d  CE11 = _mm512_set1_pd(2.518173854179933105218635e-08);
    const __m512d  CE10 = _mm512_set1_pd(2.755842049600488770111608e-07);
    const __m512d  CE9 = _mm512_set1_pd(2.755691815216689746619849e-06);
    const __m512d  CE8 = _mm512_set1_pd(2.480158383706245033920920e-05);
    const __m512d  CE7 = _mm512_set1_pd(0.0001984127043518048611841321);
    const __m512d  CE6 = _mm512_set1_pd(0.001388888889360258341755930);
    const __m512d  CE5 = _mm512_set1_pd(0.008333333332907368102819109);
    const __m512d  CE4 = _mm512_set1_pd(0.04166666666663836745814631);
    const __m512d  CE3 = _mm512_set1_pd(0.1666666666666796929434570);
    const __m512d  CE2 = _mm512_set1_pd(0.5);
    const __m512d  one = _mm512_set1_pd(1.0);
    __m512d        fexppart;
    __m512d        intpart;
    __m512d        y, p;

    x         = _mm512_max_pd(x, _mm512_div_pd(_mm512_set1_pd(std::numeric_limits<std::int32_t>::lowest()), argscale));

    y         = _mm512_mul_pd(x, argscale);

    fexppart  = ldexp(one, _mm512_cvtpd_epi32 (y));
    intpart   = _mm512_roundscale_pd(y, _MM_FROUND_TO_NEAREST_INT);

    // Extended precision arithmetics
    x         = fma(invargscale0, intpart, x);
    x         = fma(invargscale1, intpart, x);

    p         = fma(CE12, x, CE11);
    p         = fma(p, x, CE10);
    p         = fma(p, x, CE9);
    p         = fma(p, x, CE8);
    p         = fma(p, x, CE7);
    p         = fma(p, x, CE6);
    p         = fma(p, x, CE5);
    p         = fma(p, x, CE4);
    p         = fma(p, x, CE3);
    p         = fma(p, x, CE2);
    p         = fma(p, _mm512_mul_pd(x, x), x);
    x         = fma(p, fexppart, fexppart);

    return x;
}

static inline __m512d erf(__m512d x)
{
    // Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75]
    const __m512d CAP4 = _mm512_set1_pd(-0.431780540597889301512e-4);
    const __m512d CAP3 = _mm512_set1_pd(-0.00578562306260059236059);
    const __m512d CAP2 = _mm512_set1_pd(-0.028593586920219752446);
    const __m512d CAP1 = _mm512_set1_pd(-0.315924962948621698209);
    const __m512d CAP0 = _mm512_set1_pd(0.14952975608477029151);

    const __m512d CAQ5 = _mm512_set1_pd(-0.374089300177174709737e-5);
    const __m512d CAQ4 = _mm512_set1_pd(0.00015126584532155383535);
    const __m512d CAQ3 = _mm512_set1_pd(0.00536692680669480725423);
    const __m512d CAQ2 = _mm512_set1_pd(0.0668686825594046122636);
    const __m512d CAQ1 = _mm512_set1_pd(0.402604990869284362773);
    // CAQ0 == 1.0
    const __m512d CAoffset = _mm512_set1_pd(0.9788494110107421875);

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5]
    const __m512d CBP6 = _mm512_set1_pd(2.49650423685462752497647637088e-10);
    const __m512d CBP5 = _mm512_set1_pd(0.00119770193298159629350136085658);
    const __m512d CBP4 = _mm512_set1_pd(0.0164944422378370965881008942733);
    const __m512d CBP3 = _mm512_set1_pd(0.0984581468691775932063932439252);
    const __m512d CBP2 = _mm512_set1_pd(0.317364595806937763843589437418);
    const __m512d CBP1 = _mm512_set1_pd(0.554167062641455850932670067075);
    const __m512d CBP0 = _mm512_set1_pd(0.427583576155807163756925301060);
    const __m512d CBQ7 = _mm512_set1_pd(0.00212288829699830145976198384930);
    const __m512d CBQ6 = _mm512_set1_pd(0.0334810979522685300554606393425);
    const __m512d CBQ5 = _mm512_set1_pd(0.2361713785181450957579508850717);
    const __m512d CBQ4 = _mm512_set1_pd(0.955364736493055670530981883072);
    const __m512d CBQ3 = _mm512_set1_pd(2.36815675631420037315349279199);
    const __m512d CBQ2 = _mm512_set1_pd(3.55261649184083035537184223542);
    const __m512d CBQ1 = _mm512_set1_pd(2.93501136050160872574376997993);
    // CBQ0 == 1.0

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf]
    const __m512d CCP6 = _mm512_set1_pd(-2.8175401114513378771);
    const __m512d CCP5 = _mm512_set1_pd(-3.22729451764143718517);
    const __m512d CCP4 = _mm512_set1_pd(-2.5518551727311523996);
    const __m512d CCP3 = _mm512_set1_pd(-0.687717681153649930619);
    const __m512d CCP2 = _mm512_set1_pd(-0.212652252872804219852);
    const __m512d CCP1 = _mm512_set1_pd(0.0175389834052493308818);
    const __m512d CCP0 = _mm512_set1_pd(0.00628057170626964891937);

    const __m512d CCQ6 = _mm512_set1_pd(5.48409182238641741584);
    const __m512d CCQ5 = _mm512_set1_pd(13.5064170191802889145);
    const __m512d CCQ4 = _mm512_set1_pd(22.9367376522880577224);
    const __m512d CCQ3 = _mm512_set1_pd(15.930646027911794143);
    const __m512d CCQ2 = _mm512_set1_pd(11.0567237927800161565);
    const __m512d CCQ1 = _mm512_set1_pd(2.79257750980575282228);
    // CCQ0 == 1.0
    const __m512d CCoffset = _mm512_set1_pd(0.5579090118408203125);

    const __m512d one = _mm512_set1_pd(1.0);
    const __m512d two = _mm512_set1_pd(2.0);

    __m512d       xabs, x2, x4, t, t2, w, w2;
    __m512d       PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    __m512d       PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    __m512d       PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    __m512d       res_erf, res_erfcB, res_erfcC, res_erfc, res;
    __m512d       expmx2;
    __mmask8      mask, notmask_erf, mask_erf;

    // Calculate erf()
    xabs        = _mm512_abs_pd (x);
    mask_erf    = _mm512_cmp_pd_mask(xabs, one, 1);
    notmask_erf = _mm512_cmp_pd_mask(one, xabs, 2);
    x2          = _mm512_mul_pd(x, x);
    x4          = _mm512_mul_pd(x2, x2);

    PolyAP0  = fma(CAP4, x4, CAP2);
    PolyAP1  = fma(CAP3, x4, CAP1);
    PolyAP0  = fma(PolyAP0, x4, CAP0);
    PolyAP0  = fma(PolyAP1, x2, PolyAP0);

    PolyAQ1  = fma(CAQ5, x4, CAQ3);
    PolyAQ0  = fma(CAQ4, x4, CAQ2);
    PolyAQ1  = fma(PolyAQ1, x4, CAQ1);
    PolyAQ0  = fma(PolyAQ0, x4, one);
    PolyAQ0  = fma(PolyAQ1, x2, PolyAQ0);

    res_erf  = _mm512_mul_pd(PolyAP0, _mm512_maskz_rcp14_pd(mask_erf, PolyAQ0));
    res_erf  = _mm512_add_pd(CAoffset, res_erf);
    res_erf  = _mm512_mul_pd(x, res_erf);

    // Calculate erfc() in range [1,4.5]
    t       = _mm512_sub_pd(xabs, one);
    t2      = _mm512_mul_pd(t, t);

    PolyBP0  = fma(CBP6, t2, CBP4);
    PolyBP1  = fma(CBP5, t2, CBP3);
    PolyBP0  = fma(PolyBP0, t2, CBP2);
    PolyBP1  = fma(PolyBP1, t2, CBP1);
    PolyBP0  = fma(PolyBP0, t2, CBP0);
    PolyBP0  = fma(PolyBP1, t, PolyBP0);

    PolyBQ1 = fma(CBQ7, t2, CBQ5);
    PolyBQ0 = fma(CBQ6, t2, CBQ4);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ3);
    PolyBQ0 = fma(PolyBQ0, t2, CBQ2);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ1);
    PolyBQ0 = fma(PolyBQ0, t2, one);
    PolyBQ0 = fma(PolyBQ1, t, PolyBQ0);

    // The denominator polynomial can be zero outside the range
    res_erfcB  = _mm512_mul_pd(PolyBP0, _mm512_maskz_rcp14_pd(notmask_erf, PolyBQ0));

    res_erfcB = _mm512_mul_pd(res_erfcB, xabs);

    // Calculate erfc() in range [4.5,inf]
    w       = _mm512_maskz_rcp14_pd(notmask_erf, xabs);
    w2      = _mm512_mul_pd(w, w);

    PolyCP0  = fma(CCP6, w2, CCP4);
    PolyCP1  = fma(CCP5, w2, CCP3);
    PolyCP0  = fma(PolyCP0, w2, CCP2);
    PolyCP1  = fma(PolyCP1, w2, CCP1);
    PolyCP0  = fma(PolyCP0, w2, CCP0);
    PolyCP0  = fma(PolyCP1, w, PolyCP0);

    PolyCQ0  = fma(CCQ6, w2, CCQ4);
    PolyCQ1  = fma(CCQ5, w2, CCQ3);
    PolyCQ0  = fma(PolyCQ0, w2, CCQ2);
    PolyCQ1  = fma(PolyCQ1, w2, CCQ1);
    PolyCQ0  = fma(PolyCQ0, w2, one);
    PolyCQ0  = fma(PolyCQ1, w, PolyCQ0);

    expmx2   = exp( -x2 );

    // The denominator polynomial can be zero outside the range
    res_erfcC = _mm512_mul_pd(PolyCP0, _mm512_maskz_rcp14_pd(notmask_erf, PolyCQ0));
                // PolyCP0 * maskzInv(PolyCQ0, notmask_erf);
    res_erfcC = _mm512_add_pd(res_erfcC, CCoffset);
    res_erfcC = _mm512_mul_pd(res_erfcC, w);

    mask     = _mm512_cmp_pd_mask(_mm512_set1_pd(4.5), xabs, 1);
    res_erfc = _mm512_mask_blend_pd(mask, res_erfcB, res_erfcC);

    res_erfc = _mm512_mul_pd(res_erfc, expmx2);

    // erfc(x<0) = 2-erfc(|x|)
    mask     = _mm512_cmp_pd_mask(x, _mm512_setzero_pd(), 1);
    res_erfc = _mm512_mask_blend_pd(mask, res_erfc, _mm512_sub_pd(two, res_erfc));

    // Select erf() or erfc()
    res  = _mm512_mask_blend_pd(mask_erf, _mm512_sub_pd(one, res_erfc), res_erf);

    return res;
}

static inline __m512d 
erfc(__m512d x)
{
    // Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75]
    const __m512d CAP4 = _mm512_set1_pd(-0.431780540597889301512e-4);
    const __m512d CAP3 = _mm512_set1_pd(-0.00578562306260059236059);
    const __m512d CAP2 = _mm512_set1_pd(-0.028593586920219752446);
    const __m512d CAP1 = _mm512_set1_pd(-0.315924962948621698209);
    const __m512d CAP0 = _mm512_set1_pd(0.14952975608477029151);

    const __m512d CAQ5 = _mm512_set1_pd(-0.374089300177174709737e-5);
    const __m512d CAQ4 = _mm512_set1_pd(0.00015126584532155383535);
    const __m512d CAQ3 = _mm512_set1_pd(0.00536692680669480725423);
    const __m512d CAQ2 = _mm512_set1_pd(0.0668686825594046122636);
    const __m512d CAQ1 = _mm512_set1_pd(0.402604990869284362773);
    // CAQ0 == 1.0
    const __m512d CAoffset = _mm512_set1_pd(0.9788494110107421875);

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5]
    const __m512d CBP6 = _mm512_set1_pd(2.49650423685462752497647637088e-10);
    const __m512d CBP5 = _mm512_set1_pd(0.00119770193298159629350136085658);
    const __m512d CBP4 = _mm512_set1_pd(0.0164944422378370965881008942733);
    const __m512d CBP3 = _mm512_set1_pd(0.0984581468691775932063932439252);
    const __m512d CBP2 = _mm512_set1_pd(0.317364595806937763843589437418);
    const __m512d CBP1 = _mm512_set1_pd(0.554167062641455850932670067075);
    const __m512d CBP0 = _mm512_set1_pd(0.427583576155807163756925301060);
    const __m512d CBQ7 = _mm512_set1_pd(0.00212288829699830145976198384930);
    const __m512d CBQ6 = _mm512_set1_pd(0.0334810979522685300554606393425);
    const __m512d CBQ5 = _mm512_set1_pd(0.2361713785181450957579508850717);
    const __m512d CBQ4 = _mm512_set1_pd(0.955364736493055670530981883072);
    const __m512d CBQ3 = _mm512_set1_pd(2.36815675631420037315349279199);
    const __m512d CBQ2 = _mm512_set1_pd(3.55261649184083035537184223542);
    const __m512d CBQ1 = _mm512_set1_pd(2.93501136050160872574376997993);
    // CBQ0 == 1.0

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf]
    const __m512d CCP6 = _mm512_set1_pd(-2.8175401114513378771);
    const __m512d CCP5 = _mm512_set1_pd(-3.22729451764143718517);
    const __m512d CCP4 = _mm512_set1_pd(-2.5518551727311523996);
    const __m512d CCP3 = _mm512_set1_pd(-0.687717681153649930619);
    const __m512d CCP2 = _mm512_set1_pd(-0.212652252872804219852);
    const __m512d CCP1 = _mm512_set1_pd(0.0175389834052493308818);
    const __m512d CCP0 = _mm512_set1_pd(0.00628057170626964891937);

    const __m512d CCQ6 = _mm512_set1_pd(5.48409182238641741584);
    const __m512d CCQ5 = _mm512_set1_pd(13.5064170191802889145);
    const __m512d CCQ4 = _mm512_set1_pd(22.9367376522880577224);
    const __m512d CCQ3 = _mm512_set1_pd(15.930646027911794143);
    const __m512d CCQ2 = _mm512_set1_pd(11.0567237927800161565);
    const __m512d CCQ1 = _mm512_set1_pd(2.79257750980575282228);
    // CCQ0 == 1.0
    const __m512d CCoffset = _mm512_set1_pd(0.5579090118408203125);

    const __m512d one = _mm512_set1_pd(1.0);
    const __m512d two = _mm512_set1_pd(2.0);

    __m512d       xabs, x2, x4, t, t2, w, w2;
    __m512d       PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    __m512d       PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    __m512d       PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    __m512d       res_erf, res_erfcB, res_erfcC, res_erfc, res;
    __m512d       expmx2;
    __mmask8      mask, mask_erf, notmask_erf;

    // Calculate erf()
    xabs        = _mm512_abs_pd(x);
    mask_erf    = _mm512_cmp_pd_mask(xabs, one, 1); //(xabs < one);
    notmask_erf = _mm512_cmp_pd_mask(one, xabs, 2); //(one <= xabs);
    x2          = _mm512_mul_pd(x, x); //x * x;
    x4          = _mm512_mul_pd(x2, x2); //x2 * x2;

    PolyAP0  = fma(CAP4, x4, CAP2);
    PolyAP1  = fma(CAP3, x4, CAP1);
    PolyAP0  = fma(PolyAP0, x4, CAP0);
    PolyAP0  = fma(PolyAP1, x2, PolyAP0);
    PolyAQ1  = fma(CAQ5, x4, CAQ3);
    PolyAQ0  = fma(CAQ4, x4, CAQ2);
    PolyAQ1  = fma(PolyAQ1, x4, CAQ1);
    PolyAQ0  = fma(PolyAQ0, x4, one);
    PolyAQ0  = fma(PolyAQ1, x2, PolyAQ0);

    res_erf  = _mm512_mul_pd(PolyAP0, _mm512_maskz_rcp14_pd(mask_erf, PolyAQ0)); 
               //PolyAP0 * maskzInv(PolyAQ0, mask_erf);
    res_erf  = _mm512_add_pd(CAoffset, res_erf); //CAoffset + res_erf;
    res_erf  = _mm512_mul_pd(x, res_erf); //x * res_erf;

    // Calculate erfc() in range [1,4.5]
    t       = _mm512_sub_pd(xabs, one); //xabs - one;
    t2      = _mm512_mul_pd(t, t); //t * t;

    PolyBP0  = fma(CBP6, t2, CBP4);
    PolyBP1  = fma(CBP5, t2, CBP3);
    PolyBP0  = fma(PolyBP0, t2, CBP2);
    PolyBP1  = fma(PolyBP1, t2, CBP1);
    PolyBP0  = fma(PolyBP0, t2, CBP0);
    PolyBP0  = fma(PolyBP1, t, PolyBP0);

    PolyBQ1 = fma(CBQ7, t2, CBQ5);
    PolyBQ0 = fma(CBQ6, t2, CBQ4);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ3);
    PolyBQ0 = fma(PolyBQ0, t2, CBQ2);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ1);
    PolyBQ0 = fma(PolyBQ0, t2, one);
    PolyBQ0 = fma(PolyBQ1, t, PolyBQ0);

    // The denominator polynomial can be zero outside the range
    res_erfcB = _mm512_mul_pd(PolyBP0, _mm512_maskz_rcp14_pd(notmask_erf, PolyBQ0));
                //PolyBP0 * maskzInv(PolyBQ0, notmask_erf);

    res_erfcB = _mm512_mul_pd(res_erfcB, xabs); //res_erfcB * xabs;

    // Calculate erfc() in range [4.5,inf]
    w       = _mm512_maskz_rcp14_pd(_mm512_cmp_pd_mask(xabs, _mm512_setzero_pd(), 12), xabs); 
              //maskzInv(xabs, xabs != setZero());
    w2      = _mm512_mul_pd(w, w); //w * w;

    PolyCP0  = fma(CCP6, w2, CCP4);
    PolyCP1  = fma(CCP5, w2, CCP3);
    PolyCP0  = fma(PolyCP0, w2, CCP2);
    PolyCP1  = fma(PolyCP1, w2, CCP1);
    PolyCP0  = fma(PolyCP0, w2, CCP0);
    PolyCP0  = fma(PolyCP1, w, PolyCP0);

    PolyCQ0  = fma(CCQ6, w2, CCQ4);
    PolyCQ1  = fma(CCQ5, w2, CCQ3);
    PolyCQ0  = fma(PolyCQ0, w2, CCQ2);
    PolyCQ1  = fma(PolyCQ1, w2, CCQ1);
    PolyCQ0  = fma(PolyCQ0, w2, one);
    PolyCQ0  = fma(PolyCQ1, w, PolyCQ0);

    expmx2   = exp( -x2 );

    // The denominator polynomial can be zero outside the range
    res_erfcC = _mm512_mul_pd(PolyCP0, _mm512_maskz_rcp14_pd(notmask_erf, PolyCQ0));
                //PolyCP0 * maskzInv(PolyCQ0, notmask_erf);
    res_erfcC = _mm512_add_pd(res_erfcC, CCoffset); 
                //res_erfcC + CCoffset;
    res_erfcC = _mm512_mul_pd(res_erfcC, w);
                //res_erfcC * w;

    mask     = _mm512_cmp_pd_mask(_mm512_set1_pd(4.5), xabs, 1); //(__m512d(4.5) < xabs);
    res_erfc = _mm512_mask_blend_pd(mask, res_erfcB, res_erfcC);//blend(res_erfcB, res_erfcC, mask);

    res_erfc = _mm512_mul_pd(res_erfc, expmx2);//res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    mask     = _mm512_cmp_pd_mask(x, _mm512_setzero_pd(), 1);//(x < setZero());
    res_erfc = _mm512_mask_blend_pd(mask, res_erfc, _mm512_sub_pd(two, res_erfc));//blend(res_erfc, two - res_erfc, mask);

    // Select erf() or erfc()
    res  = _mm512_mask_blend_pd(mask_erf, res_erfc, _mm512_sub_pd(one, res_erf));
           //blend(res_erfc, one - res_erf, mask_erf);

    return res;
}

#endif

#endif