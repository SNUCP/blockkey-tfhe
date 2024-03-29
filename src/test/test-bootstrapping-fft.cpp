#include "lwekey.h"
#include "lweparams.h"
#include "lwesamples.h"
#include "polynomials.h"
#include "tfhe.h"
#include "tgsw.h"
#include "tlwe.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <sys/time.h>

using namespace std;

// **********************************************************************************
// ********************************* MAIN
// *******************************************
// **********************************************************************************

void dieDramatically(string message) {
  cerr << message << endl;
  abort();
}

void fillRandom(LweSample *result, const LweParams *params) {
  const int32_t n = params->n;
  for (int32_t i = 0; i < n; i++)
    result->a[i] = uniformTorus32_distrib(generator);
  result->b = uniformTorus32_distrib(generator);
  result->current_variance = 0.2;
}

// EXPORT void tLweExtractKey(LweKey* result, const TLweKey* key); //TODO:
// change the name and put in a .h EXPORT void
// tfhe_createLweBootstrappingKeyFFT(LweBootstrappingKeyFFT* bk, const LweKey*
// key_in, const TGswKey* rgsw_key); EXPORT void tfhe_bootstrapFFT(LweSample*
// result, const LweBootstrappingKeyFFT* bk, Torus32 mu1, Torus32 mu0, const
// LweSample* x);

#ifndef NDEBUG
extern const TLweKey *debug_accum_key;
extern const LweKey *debug_extract_key;
extern const LweKey *debug_in_key;
#endif

int32_t main(int32_t argc, char **argv) {
#ifndef NDEBUG
  cout << "DEBUG MODE!" << endl;
#endif
  const int32_t nb_samples = 100;
  const Torus32 mu_boot = modSwitchToTorus32(1, 8);

  // generate params
  int32_t minimum_lambda = 100;
  TFheGateBootstrappingParameterSet *params =
      new_default_gate_bootstrapping_parameters(minimum_lambda);
  const LweParams *in_out_params = params->in_out_params;
  // generate the secret keyset
  TFheGateBootstrappingSecretKeySet *keyset =
      new_random_gate_bootstrapping_secret_keyset(params);

  // generate input samples
  LweSample *test_in = new_LweSample_array(nb_samples, in_out_params);
  for (int32_t i = 0; i < nb_samples; ++i) {
    lweSymEncrypt(test_in + i, modSwitchToTorus32(i, nb_samples), 0.01,
                  keyset->lwe_key);
  }
  // output samples
  LweSample *test_out = new_LweSample_array(nb_samples, in_out_params);

  clock_t begin, end;

  // bootstrap input samples
  cout << "starting bootstrapping..." << endl;
  begin = clock();
  for (int32_t i = 0; i < nb_samples; ++i) {
    tfhe_bootstrap_FFT(test_out + i, keyset->cloud.bkFFT, mu_boot, test_in + i);
  }
  end = clock();
  cout << "finished " << nb_samples << " bootstrappings" << endl;
  cout << "time per bootstrapping (microsecs)... "
       << (end - begin) / double(nb_samples) << endl;

  delete_LweSample_array(nb_samples, test_in);

  /** keyswitch **/
  test_in = new_LweSample_array(
      nb_samples, &keyset->cloud.bkFFT->accum_params->extracted_lweparams);

  for (int32_t i = 0; i < nb_samples; ++i) {
    fillRandom(test_in + i,
               &keyset->cloud.bkFFT->accum_params->extracted_lweparams);
  }

  cout << "starting key-switching..." << endl;
  begin = clock();
  for (int32_t i = 0; i < nb_samples; ++i) {
    lweKeySwitch(test_out + i, keyset->cloud.bkFFT->ks, test_in + i);
  }
  end = clock();
  cout << "finished " << nb_samples << " key-switching" << endl;
  cout << "time per key-switching (microsecs)... "
       << (end - begin) / double(nb_samples) << endl;

  delete_LweSample_array(nb_samples, test_out);
  delete_LweSample_array(nb_samples, test_in);

  delete_gate_bootstrapping_secret_keyset(keyset);
  delete_gate_bootstrapping_parameters(params);

  /** sparse boot **/
  params = new_sparse_gate_bootstrapping_parameters();
  in_out_params = params->in_out_params;
  // generate the secret keyset
  keyset = new_random_sparse_bootstrapping_secret_keyset(params);

  // generate input samples
  test_in = new_LweSample_array(nb_samples, in_out_params);
  for (int32_t i = 0; i < nb_samples; ++i) {
    lweSymEncrypt(test_in + i, modSwitchToTorus32(i, nb_samples), 0.01,
                  keyset->lwe_key);
  }
  // output samples
  test_out = new_LweSample_array(nb_samples, in_out_params);

  // bootstrap input samples
  cout << "starting sparse bootstrapping..." << endl;
  begin = clock();
  for (int32_t i = 0; i < nb_samples; ++i) {
    tfhe_sparseBootstrap_FFT(test_out + i, keyset->params->hw,
                             keyset->cloud.bkFFT, mu_boot, test_in + i);
  }

  end = clock();
  cout << "finished " << nb_samples << " sparse bootstrappings" << endl;
  cout << "time per sparse bootstrapping (microsecs)... "
       << (end - begin) / double(nb_samples) << endl;

  delete_LweSample_array(nb_samples, test_in);

  /** sparse keyswitch **/
  test_in = new_LweSample_array(
      nb_samples, &keyset->cloud.bkFFT->accum_params->extracted_lweparams);

  for (int32_t i = 0; i < nb_samples; ++i) {
    fillRandom(test_in + i,
               &keyset->cloud.bkFFT->accum_params->extracted_lweparams);
  }

  cout << "starting sparse key-switching..." << endl;
  begin = clock();
  for (int32_t i = 0; i < nb_samples; ++i) {
    lweSparseKeySwitch(test_out + i, keyset->cloud.bkFFT->ks, test_in + i);
  }
  end = clock();
  cout << "finished " << nb_samples << " sparse key-switching" << endl;
  cout << "time per sparse key-switching (microsecs)... "
       << (end - begin) / double(nb_samples) << endl;

  delete_LweSample_array(nb_samples, test_out);
  delete_LweSample_array(nb_samples, test_in);

  delete_gate_bootstrapping_secret_keyset(keyset);
  delete_gate_bootstrapping_parameters(params);

  return 0;
}
