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
#ifndef PITCH_MATRIX_CBLAS_WRAPPERS_H_
#define PITCH_MATRIX_CBLAS_WRAPPERS_H_ 1


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


/**This function is a matrix-level scale operation. In brief, 
 * the expression is "A = alpha * A(Matrix level)".
 *
 * @param [in,out] *data  It always be a MatrixBase pointer in kaldi. It is the
 *                        entity of data.
 * @param [in] alpha  Scaling factor for the values in "data".
 * @param [in] num_rows  Number of rows in matrix "data".   
 * @param [in] num_cols  Number of columns in matrix "data".
 * @param [in] stride  The distance in memory betwwen each row.
 */
template<typename Real>
inline void blas_scale(Real* data, const Real alpha, const MatrixIndexT num_rows, 
               const MatrixIndexT num_cols, const MatrixIndexT stride) {
  if (alpha == 1.0) { 
    return;
  } else { //alpha != 1.0
    for(MatrixIndexT r = 0; r < num_rows; r++) {
      for(MatrixIndexT c = 0; c < num_cols; c++) {
        data[r*stride + c] *= alpha;
      }  
    }
  }
}


// template instantiation(float and double)
template
void blas_scale(float* data, const float alpha, const MatrixIndexT num_rows,
                       const MatrixIndexT num_cols, const MatrixIndexT stride);
template
void blas_scale(double* data, const double alpha, const MatrixIndexT num_rows,
                       const MatrixIndexT num_cols, const MatrixIndexT stride);


/**This function does a constant times a vector plus a vector. In brief,
 * the expression is ori = ori + alpha * additional(vector level).
 *
 * @param [in] num  It is INTEGER, number of elements need to be dealed.
 * @param [in] alpha  Scaling factor for the values in "additional_data".
 * @param [in] *additional_data  It always be a VectorBase pointer. It contains
 *                               the entity of data which will multiply by a 
 *                               constant.
 * @param [in] additional_inc  Stride with "additional_data". Each 
 *                             "additional_inc"th element will be used.
 * @param [in,out] *ori_data  The entity of data. The values will be used later.
 *                            It add the result of multiplication.
 * @param [in] ori_inc  Stride with "ori_data". Each "ori_inc"th element will
 *                      be used.
 */
template<typename Real>
inline void blas_axpy(const MatrixIndexT num, const Real alpha,
                      const Real* additional_data, const MatrixIndexT additional_inc,
                      Real* ori_data, const MatrixIndexT ori_inc) {
  for(MatrixIndexT i = 0; i < num; i++) {
    if (alpha != 1.0) {
      ori_data[i * ori_inc] += alpha * additional_data[i * additional_inc];
    } else {
      ori_data[i * ori_inc] += additional_data[i * additional_inc];
    }
  }
}


// template instantiation(float and double)
template
void blas_axpy(const MatrixIndexT num, const float alpha, const float* additional_data,
               const MatrixIndexT additional_inc, float* ori_data, const MatrixIndexT ori_inc);
template
void blas_axpy(const MatrixIndexT num, const double alpha, const double* additional_data,
               const MatrixIndexT additional_inc, double* ori_data, const MatrixIndexT ori_inc);


/**The function does inner product of two vectors.
 * 
 * @param [in] num  The number of elements need to be dealed.
 * @param [in] X  Input vector X.
 * @param [in] incX  Stride within X. Each "incX"th element is used.
 * @param [in] Y  Input vector Y.
 * @param [in] incY  Stride within Y. Each "incY"th element is used.
 *
 * @return Return the result of inner product.
 */
template<typename Real>
inline Real blas_dot(const MatrixIndexT num, const Real *const X, const int incX,
                     const Real *const Y, const MatrixIndexT incY){
  Real result = 0.0;
  for(MatrixIndexT i = 0; i < num; i++) {
    result = result +  X[i * incX] * Y[i * incY];
  }
  return result;
}


// template instantiation(float and double)
template
float blas_dot(const MatrixIndexT num, const float *const X, const MatrixIndexT incX,
               const float *const Y, const MatrixIndexT incY);
template
double blas_dot(const MatrixIndexT num, const double *const X, const MatrixIndexT incX,
                const double *const Y, const MatrixIndexT incY);


/**The function does vector-level scale operation. In brief, the expression is
 * (*data) = alpha * (*data).[vector level]
 *
 * @param [in] num  The number of elements need to be operated.
 * @param [in] alpha  Scaling factor for the values in "data".
 * @param [in,out] *data  It always be a VectorBase pointer in kaldi. It contains
 *                        the entity of data.
 * @param [in] inc  Stride within "data". Each "inc"th element will be used.
 */
template<typename Real>
inline void blas_scale(const MatrixIndexT num, const float alpha, Real *data, const int inc){
  for(MatrixIndexT i = 0; i < num; i++) {
    if (alpha != 1.0) data[i * inc] *= alpha;
  }
}


// template instantiation(float and double)
template
void blas_scale(const MatrixIndexT num, const float alpha, float *data, const int inc);
template
void blas_scale(const MatrixIndexT num, const float alpha, double *data, const int inc);
}
// namespace kaldi

#endif
