/*
 * Bootstrapping FFT functions
 */

#ifndef TFHE_TEST_ENVIRONMENT

#include "tfhe.h"
#include <cassert>
#include <iostream>

using namespace std;
#define INCLUDE_ALL
#else
#undef EXPORT
#define EXPORT
#endif

void tfhe_sparseMuxRotate_FFT(TLweSample *result, const TLweSample *accum,
                              const TGswSampleFFT *bki, const int32_t *bara,
                              const int32_t d, const TGswParams *bk_params) {
  // ACC = BKi*[(X^barai-1)*ACC]+ACC
  // temp = (X^barai-1)*ACC
  tLweCopy(result, accum, bk_params->tlwe_params);
  // temp *= BKi

  tGswFFTExternMulToTLweHoisting(result, bki, bara, d, bk_params);

  // ACC += temp
  tLweAddTo(result, accum, bk_params->tlwe_params);
}

#if defined INCLUDE_ALL || defined INCLUDE_TFHE_BLIND_ROTATE_FFT
#undef INCLUDE_TFHE_BLIND_ROTATE_FFT
/**
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
EXPORT void tfhe_sparseBlindRotate_FFT(TLweSample *accum,
                                       const TGswSampleFFT *bkFFT,
                                       const int32_t *bara, const int32_t n,
                                       const int32_t hw,
                                       const TLweParams *params,
                                       const TGswParams *bk_params) {

  const int32_t d = n / hw;
  TLweSample *temp = new_TLweSample(bk_params->tlwe_params);
  TLweSample *temp2 = temp;
  TLweSample *temp3 = accum;

  for (int32_t i = 0; i < hw; i++) {

    tfhe_sparseMuxRotate_FFT(temp2, temp3, bkFFT + i * d, bara + i * d, d,
                             bk_params);
    swap(temp2, temp3);
  }
  if (temp3 != accum) {
    tLweCopy(accum, temp3, bk_params->tlwe_params);
  }

  delete_TLweSample(temp);
}
#endif

#if defined INCLUDE_ALL || defined INCLUDE_TFHE_BLIND_ROTATE_AND_EXTRACT_FFT
#undef INCLUDE_TFHE_BLIND_ROTATE_AND_EXTRACT_FFT
/**
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
EXPORT void tfhe_sparseBlindRotateAndExtract_FFT(
    LweSample *result, const TorusPolynomial *v, const TGswSampleFFT *bk,
    const int32_t barb, const int32_t *bara, const int32_t n, const int32_t hw,
    const TGswParams *bk_params) {

  const TLweParams *accum_params = bk_params->tlwe_params;
  const LweParams *extract_params = &accum_params->extracted_lweparams;
  const int32_t N = accum_params->N;
  const int32_t _2N = 2 * N;

  // Test polynomial
  TorusPolynomial *testvectbis = new_TorusPolynomial(N);
  // Accumulator
  TLweSample *acc = new_TLweSample(accum_params);

  int32_t temp = (_2N - barb) % _2N;
  /*
  const int32_t d = n / hw;

  for (int32_t i = 0; i < hw; i++) {
    temp += bara[i * d + d - 1];
  }
  */

  // testvector = X^{temp}*v
  if (temp != 0)
    torusPolynomialMulByXai(testvectbis, temp, v);
  else
    torusPolynomialCopy(testvectbis, v);

  tLweNoiselessTrivial(acc, testvectbis, accum_params);
  // Blind rotation
  tfhe_sparseBlindRotate_FFT(acc, bk, bara, n, hw, accum_params, bk_params);
  // Extraction
  tLweExtractLweSample(result, acc, extract_params, accum_params);

  delete_TLweSample(acc);
  delete_TorusPolynomial(testvectbis);
}
#endif

#if defined INCLUDE_ALL || defined INCLUDE_TFHE_BOOTSTRAP_WO_KS_FFT
#undef INCLUDE_TFHE_BOOTSTRAP_WO_KS_FFT
/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
EXPORT void tfhe_sparseBootstrap_woKS_FFT(LweSample *result, const int32_t hw,
                                          const LweBootstrappingKeyFFT *bk,
                                          Torus32 mu, const LweSample *x) {

  const TGswParams *bk_params = bk->bk_params;
  const TLweParams *accum_params = bk->accum_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t N = accum_params->N;
  const int32_t Nx2 = 2 * N;
  const int32_t n = in_params->n;

  TorusPolynomial *testvect = new_TorusPolynomial(N);
  int32_t *bara = new int32_t[N];

  // Modulus switching
  int32_t barb = modSwitchFromTorus32(x->b, Nx2);
  for (int32_t i = 0; i < n; i++) {
    bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
  }

  // the initial testvec = [mu,mu,mu,...,mu]
  for (int32_t i = 0; i < N; i++)
    testvect->coefsT[i] = mu;

  // Bootstrapping rotation and extraction
  tfhe_sparseBlindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara,
                                       n, hw, bk_params);

  delete[] bara;
  delete_TorusPolynomial(testvect);
}
#endif

#if defined INCLUDE_ALL || defined INCLUDE_TFHE_BOOTSTRAP_FFT
#undef INCLUDE_TFHE_BOOTSTRAP_FFT
/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
EXPORT void tfhe_sparseBootstrap_FFT(LweSample *result, const int32_t hw,
                                     const LweBootstrappingKeyFFT *bk,
                                     Torus32 mu, const LweSample *x) {

  LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

  tfhe_sparseBootstrap_woKS_FFT(u, hw, bk, mu, x);
  // Key switching
  lweSparseKeySwitch(result, bk->ks, u);

  delete_LweSample(u);
}
#endif
