// matrix/cblas-wrappers.h

// Copyright 2017  Johns Hopkins University (author: Daniel Povey);
//                 Hang Lyu

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//  http://www.apache.org/licenses/LICENSE-2.0

// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.
#ifndef KALDI_MATRIX_CBLAS_WRAPPERS_H_
#define KALDI_MATRIX_CBLAS_WRAPPERS_H_ 1


#include <limits>
#include "matrix/kaldi-vector.h"
#include "matrix/kaldi-matrix.h"
#include "matrix/matrix-functions.h"

// Do not include this file directly.  It is to be included
// by .cc files in this directory.

namespace kaldi {

inline void mul_elements(
    const MatrixIndexT dim,
    const double *a,
    double *b) { // does b *= a, elementwise.
  double c1, c2, c3, c4;
  MatrixIndexT i;
  for (i = 0; i + 4 <= dim; i += 4) {
    c1 = a[i] * b[i];
    c2 = a[i+1] * b[i+1];
    c3 = a[i+2] * b[i+2];
    c4 = a[i+3] * b[i+3];
    b[i] = c1;
    b[i+1] = c2;
    b[i+2] = c3;
    b[i+3] = c4;
  }
  for (; i < dim; i++)
    b[i] *= a[i];
}

inline void mul_elements(
    const MatrixIndexT dim,
    const float *a,
    float *b) { // does b *= a, elementwise.
  float c1, c2, c3, c4;
  MatrixIndexT i;
  for (i = 0; i + 4 <= dim; i += 4) {
    c1 = a[i] * b[i];
    c2 = a[i+1] * b[i+1];
    c3 = a[i+2] * b[i+2];
    c4 = a[i+3] * b[i+3];
    b[i] = c1;
    b[i+1] = c2;
    b[i+2] = c3;
    b[i+3] = c4;
  }
  for (; i < dim; i++)
    b[i] *= a[i];
}

// A = alpha * A(Matrix level)
template<typename Real>
inline void blas_scale(Real* data, const Real alpha, const MatrixIndexT num_rows, 
               const MatrixIndexT num_cols, const MatrixIndexT stride) {
  for(MatrixIndexT r = 0; r < num_rows; r++) {
    for(MatrixIndexT c = 0; c < num_cols; c++) {
      data[r * stride + c] *= alpha;
    }  
  }
}

template
void blas_scale(float* data, const float alpha, const MatrixIndexT num_rows,
                       const MatrixIndexT num_cols, const MatrixIndexT stride);
template
void blas_scale(double* data, const double alpha, const MatrixIndexT num_rows,
                       const MatrixIndexT num_cols, const MatrixIndexT stride);

// ori = ori + alpha * additional(vector level)
template<typename Real>
inline void blas_axpy(const MatrixIndexT num, const Real alpha,
                      const Real* additional_data, const MatrixIndexT additional_inc,
                      Real* ori_data, const MatrixIndexT ori_inc) {
  for(MatrixIndexT i = 0; i < num; i++) {
    ori_data[i * ori_inc] += alpha * additional_data[i * additional_inc];
  }
}

template
void blas_axpy(const MatrixIndexT num, const float alpha, const float* additional_data,
               const MatrixIndexT additional_inc, float* ori_data, const MatrixIndexT ori_inc);
template
void blas_axpy(const MatrixIndexT num, const double alpha, const double* additional_data,
               const MatrixIndexT additional_inc, double* ori_data, const MatrixIndexT ori_inc);

template<typename Real>
inline Real blas_dot(const MatrixIndexT N, const Real *const X, const int incX,
                     const Real *const Y, const MatrixIndexT incY){
  Real result = 0;
  for(MatrixIndexT i = 0; i < N; i++) {
    result = result +  X[i * incX] * Y[i * incY];
  }
  return result;
}

template
float blas_dot(const MatrixIndexT N, const float *const X, const MatrixIndexT incX,
               const float *const Y, const MatrixIndexT incY);
template
double blas_dot(const MatrixIndexT N, const double *const X, const MatrixIndexT incX,
                const double *const Y, const MatrixIndexT incY);

// *data = alpha * (*data)(vector level)
template<typename Real>
inline void blas_scale(const int N, const float alpha, Real *data, const int inc){
  for(MatrixIndexT i = 0; i < N; i++) {
    data[i * inc] *= alpha;
  }
}

template
void blas_scale(const int N, const float alpha, float *data, const int inc);
template
void blas_scale(const int N, const float alpha, double *data, const int inc);

}
// namespace kaldi

#endif
