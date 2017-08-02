// matrix/kaldi-vector.cc

// Copyright 2009-2011  Microsoft Corporation;  Lukas Burget;
//                      Saarland University;   Go Vivace Inc.;  Ariya Rastrow;
//                      Petr Schwarz;  Yanmin Qian;  Jan Silovsky;
//                      Haihua Xu; Wei Shi
//                2015  Guoguo Chen


// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <string>
#include "matrix/kaldi-vector.h"
#include "matrix/cblas-wrappers.h"

namespace kaldi {

template<typename Real>
void VectorBase<Real>::SetZero() {
  std::memset(data_, 0, dim_ * sizeof(Real));
}

template<typename Real>
bool VectorBase<Real>::IsZero(Real cutoff) const {
  Real abs_max = 0.0;
  for (MatrixIndexT i = 0; i < Dim(); i++)
    abs_max = std::max(std::abs(data_[i]), abs_max);
  return (abs_max <= cutoff);
}

template<typename Real>
void VectorBase<Real>::Set(Real f) {
  // Why not use memset here?
  for (MatrixIndexT i = 0; i < dim_; i++) { data_[i] = f; }
}

template<typename Real>
void VectorBase<Real>::SetRandn() {
  kaldi::RandomState rstate;
  MatrixIndexT last = (Dim() % 2 == 1) ? Dim() - 1 : Dim();
  for (MatrixIndexT i = 0; i < last; i += 2) {
    kaldi::RandGauss2(data_ + i, data_ + i +1, &rstate);
  }
  if (Dim() != last) data_[last] = static_cast<Real>(kaldi::RandGauss(&rstate));
}

/// Copy data from another vector
template<typename Real>
void VectorBase<Real>::CopyFromVec(const VectorBase<Real> &v) {
  KALDI_ASSERT(Dim() == v.Dim());
  if (data_ != v.data_) {
    std::memcpy(this->data_, v.data_, dim_ * sizeof(Real));
  }
}

template<typename Real>
template<typename OtherReal>
void VectorBase<Real>::CopyFromVec(const VectorBase<OtherReal> &other) {
  KALDI_ASSERT(dim_ == other.Dim());
  Real * __restrict__  ptr = data_;
  const OtherReal * __restrict__ other_ptr = other.Data();
  for (MatrixIndexT i = 0; i < dim_; i++)
    ptr[i] = other_ptr[i];
}

template void VectorBase<float>::CopyFromVec(const VectorBase<double> &other);
template void VectorBase<double>::CopyFromVec(const VectorBase<float> &other);

template<typename Real>
MatrixIndexT VectorBase<Real>::ApplyFloor(Real floor_val) {
  MatrixIndexT num_floored = 0;
  for (MatrixIndexT i = 0; i < dim_; i++) {
    if (data_[i] < floor_val) {
      data_[i] = floor_val;
      num_floored++;
    }
  }
  return num_floored;
}

template<typename Real>
MatrixIndexT VectorBase<Real>::ApplyFloor(const VectorBase<Real> &floor_vec) {
  KALDI_ASSERT(floor_vec.Dim() == dim_);
  MatrixIndexT num_floored = 0;
  for (MatrixIndexT i = 0; i < dim_; i++) {
    if (data_[i] < floor_vec(i)) {
      data_[i] = floor_vec(i);
      num_floored++;
    }
  }
  return num_floored;
}

// takes elements to a power.  Throws exception if could not (but only for power != 1 and power != 2).
template<typename Real>
void VectorBase<Real>::ApplyPow(Real power) {
  if (power == 1.0) return;
  if (power == 2.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      data_[i] = data_[i] * data_[i];
  } else if (power == 0.5) {
    for (MatrixIndexT i = 0; i < dim_; i++) {
      if (!(data_[i] >= 0.0))
        KALDI_ERR << "Cannot take square root of negative value "
                  << data_[i];
      data_[i] = std::sqrt(data_[i]);
    }
  } else {
    for (MatrixIndexT i = 0; i < dim_; i++) {
      data_[i] = pow(data_[i], power);
      if (data_[i] == HUGE_VAL) {  // HUGE_VAL is what errno returns on error.
        KALDI_ERR << "Could not raise element "  << i << " to power "
                  << power << ": returned value = " << data_[i];
      }
    }
  }
}

// Computes the p-th norm. Throws exception if could not.
template<typename Real>
Real VectorBase<Real>::Norm(Real p) const {
  KALDI_ASSERT(p >= 0.0);
  Real sum = 0.0;
  if (p == 0.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      if (data_[i] != 0.0) sum += 1.0;
    return sum;
  } else if (p == 1.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      sum += std::abs(data_[i]);
    return sum;
  } else if (p == 2.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      sum += data_[i] * data_[i];
    return std::sqrt(sum);
  } else if (p == std::numeric_limits<Real>::infinity()){
    for (MatrixIndexT i = 0; i < dim_; i++)
      sum = std::max(sum, std::abs(data_[i]));
    return sum;
  } else {
    Real tmp;
    bool ok = true;
    for (MatrixIndexT i = 0; i < dim_; i++) {
      tmp = pow(std::abs(data_[i]), p);
      if (tmp == HUGE_VAL) // HUGE_VAL is what pow returns on error.
        ok = false;
      sum += tmp;
    }
    tmp = pow(sum, static_cast<Real>(1.0/p));
    KALDI_ASSERT(tmp != HUGE_VAL); // should not happen here.
    if (ok) {
      return tmp;
    } else {
      Real maximum = this->Max(), minimum = this->Min(),
          max_abs = std::max(maximum, -minimum);
      KALDI_ASSERT(max_abs > 0); // Or should not have reached here.
      Vector<Real> tmp(*this);
      tmp.Scale(1.0 / max_abs);
      return tmp.Norm(p) * max_abs;
    }
  }
}

template<typename Real>
bool VectorBase<Real>::ApproxEqual(const VectorBase<Real> &other, float tol) const {
  if (dim_ != other.dim_) KALDI_ERR << "ApproxEqual: size mismatch "
                                    << dim_ << " vs. " << other.dim_;
  KALDI_ASSERT(tol >= 0.0);
  if (tol != 0.0) {
    Vector<Real> tmp(*this);
    tmp.AddVec(-1.0, other);
    return (tmp.Norm(2.0) <= static_cast<Real>(tol) * this->Norm(2.0));
  } else { // Test for exact equality.
    const Real *data = data_;
    const Real *other_data = other.data_;
    for (MatrixIndexT dim = dim_, i = 0; i < dim; i++)
      if (data[i] != other_data[i]) return false;
    return true;
  }
}

template<typename Real>
void VectorBase<Real>::InvertElements() {
  for (MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] = static_cast<Real>(1 / data_[i]);
  }
}

template<typename Real>
template<typename OtherReal>
void VectorBase<Real>::AddVec(const Real alpha, const VectorBase<OtherReal> &v) {
  KALDI_ASSERT(dim_ == v.dim_);
  // remove __restrict__ if it causes compilation problems.
  Real *__restrict__ data = data_;
  OtherReal *__restrict__ other_data = v.data_;
  MatrixIndexT dim = dim_;
  if (alpha != 1.0)
    for (MatrixIndexT i = 0; i < dim; i++)
      data[i] += alpha * other_data[i];
  else
    for (MatrixIndexT i = 0; i < dim; i++)
      data[i] += other_data[i];
}

template
void VectorBase<float>::AddVec(const float alpha, const VectorBase<double> &v);
template
void VectorBase<double>::AddVec(const double alpha, const VectorBase<float> &v);

template<typename Real>
void VectorBase<Real>::AddVec2(const Real alpha, const VectorBase<Real> &v) {
  KALDI_ASSERT(dim_ == v.dim_);
  for (MatrixIndexT i = 0; i < dim_; i++)
    data_[i] += alpha * v.data_[i] * v.data_[i];
}

template<typename Real>
void VectorBase<Real>::AddMatVec(const Real alpha,
                                  const MatrixBase<Real> &M,
                                  MatrixTransposeType trans,
                                  const VectorBase<Real> &v,
                                  const Real beta) {
  KALDI_ASSERT((trans == kNoTrans && M.NumCols() == v.dim_ && M.NumRows() == dim_)
               || (trans == kTrans && M.NumRows() == v.dim_ && M.NumCols() == dim_));
  KALDI_ASSERT(&v != this);
  //cblas_Xgemv(trans, M.NumRows(), M.NumCols(), alpha, M.Data(), M.Stride(),
  //            v.Data(), 1, beta, data_, 1);
  if (trans == kNoTrans) {
    MatrixIndexT dim = this->dim_, mat_stride = M.Stride();
    Real *data = this->data_;
    const Real *mat_data = M.Data(), *v_data = v.Data();
    for (MatrixIndexT i = 0; i < dim; i++, mat_data += mat_stride, data++) {
      *data = beta * (*data) + alpha * blas_dot(dim, mat_data, 1, v_data, 1);
    }
  } else { //trans == kTrans
    MatrixIndexT dim = this->dim_, mat_stride = M.Stride();
    Real *data = this->data_;
    const Real *mat_data = M.Data(), *v_data = v.Data();
    for (MatrixIndexT i = 0; i < dim; i++, mat_data++, data++) {
      *data = beta * (*data) + alpha * blas_dot(dim, mat_data, mat_stride,
                                                v_data, 1);
    }
  }
}

template<typename Real>
void VectorBase<Real>::MulElements(const VectorBase<Real> &v) {
  KALDI_ASSERT(dim_ == v.dim_);
  for (MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] *= v.data_[i];
  }
}

template<typename Real>
template<typename OtherReal>
void VectorBase<Real>::MulElements(const VectorBase<OtherReal> &v) {
  KALDI_ASSERT(dim_ == v.Dim());
  const OtherReal *other_ptr = v.Data();
  for (MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] *= other_ptr[i];
  }
}
// instantiate template.
template
void VectorBase<float>::MulElements(const VectorBase<double> &v);
template
void VectorBase<double>::MulElements(const VectorBase<float> &v);

template<typename Real>
void VectorBase<Real>::Add(Real c) {
  for (MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] += c;
  }
}

template<typename Real>
void VectorBase<Real>::AddVecVec(Real alpha, const VectorBase<Real> &v,
                                 const VectorBase<Real> &r, Real beta) {
  KALDI_ASSERT(v.data_ != this->data_ && r.data_ != this->data_);
  // We pretend that v is a band-diagonal matrix.
  KALDI_ASSERT(dim_ == v.dim_ && dim_ == r.dim_);
  // (*this) = (beta)(*this) + alpha * v * r
  for(MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] = beta * data_[i] + alpha * blas_dot(dim_, v.data_, 1, r.data_, 1);
  } 
}

template<typename Real>
void VectorBase<Real>::Scale(Real alpha) {
  blas_scale(dim_, alpha, data_, 1);
}

template<typename Real>
Real VectorBase<Real>::Max() const {
  Real ans = - std::numeric_limits<Real>::infinity();
  const Real *data = data_;
  MatrixIndexT i, dim = dim_;
  for (i = 0; i + 4 <= dim; i += 4) {
    Real a1 = data[i], a2 = data[i+1], a3 = data[i+2], a4 = data[i+3];
    if (a1 > ans || a2 > ans || a3 > ans || a4 > ans) {
      Real b1 = (a1 > a2 ? a1 : a2), b2 = (a3 > a4 ? a3 : a4);
      if (b1 > ans) ans = b1;
      if (b2 > ans) ans = b2;
    }
  }
  for (; i < dim; i++)
    if (data[i] > ans) ans = data[i];
  return ans;
}

template<typename Real>
Real VectorBase<Real>::Max(MatrixIndexT *index_out) const {
  if (dim_ == 0) KALDI_ERR << "Empty vector";
  Real ans = - std::numeric_limits<Real>::infinity();
  MatrixIndexT index = 0;
  const Real *data = data_;
  MatrixIndexT i, dim = dim_;
  for (i = 0; i + 4 <= dim; i += 4) {
    Real a1 = data[i], a2 = data[i+1], a3 = data[i+2], a4 = data[i+3];
    if (a1 > ans || a2 > ans || a3 > ans || a4 > ans) {
      if (a1 > ans) { ans = a1; index = i; }
      if (a2 > ans) { ans = a2; index = i + 1; }
      if (a3 > ans) { ans = a3; index = i + 2; }
      if (a4 > ans) { ans = a4; index = i + 3; }
    }
  }
  for (; i < dim; i++)
    if (data[i] > ans) { ans = data[i]; index = i; }
  *index_out = index;
  return ans;
}

template<typename Real>
Real VectorBase<Real>::Min() const {
  Real ans = std::numeric_limits<Real>::infinity();
  const Real *data = data_;
  MatrixIndexT i, dim = dim_;
  for (i = 0; i + 4 <= dim; i += 4) {
    Real a1 = data[i], a2 = data[i+1], a3 = data[i+2], a4 = data[i+3];
    if (a1 < ans || a2 < ans || a3 < ans || a4 < ans) {
      Real b1 = (a1 < a2 ? a1 : a2), b2 = (a3 < a4 ? a3 : a4);
      if (b1 < ans) ans = b1;
      if (b2 < ans) ans = b2;
    }
  }
  for (; i < dim; i++)
    if (data[i] < ans) ans = data[i];
  return ans;
}

template<typename Real>
Real VectorBase<Real>::Min(MatrixIndexT *index_out) const {
  if (dim_ == 0) KALDI_ERR << "Empty vector";
  Real ans = std::numeric_limits<Real>::infinity();
  MatrixIndexT index = 0;
  const Real *data = data_;
  MatrixIndexT i, dim = dim_;
  for (i = 0; i + 4 <= dim; i += 4) {
    Real a1 = data[i], a2 = data[i+1], a3 = data[i+2], a4 = data[i+3];
    if (a1 < ans || a2 < ans || a3 < ans || a4 < ans) {
      if (a1 < ans) { ans = a1; index = i; }
      if (a2 < ans) { ans = a2; index = i + 1; }
      if (a3 < ans) { ans = a3; index = i + 2; }
      if (a4 < ans) { ans = a4; index = i + 3; }
    }
  }
  for (; i < dim; i++)
    if (data[i] < ans) { ans = data[i]; index = i; }
  *index_out = index;
  return ans;
}

template<typename Real>
Real VectorBase<Real>::Sum() const {
  double sum = 0.0;
  for (MatrixIndexT i = 0; i < dim_; i++) { sum += data_[i]; }
  return sum;
}

template<typename Real>
void VectorBase<Real>::AddRowSumMat(Real alpha, const MatrixBase<Real> &M, Real beta) {
  KALDI_ASSERT(dim_ == M.NumCols());
  MatrixIndexT num_rows = M.NumRows(), stride = M.Stride(), dim = dim_;
  Real *data = data_;

  // implement the function according to a dimension cutoff for computation efficiency
  if (num_rows <= 64) {
    blas_scale(dim, beta, data, 1);
    const Real *m_data = M.Data();
    for (MatrixIndexT i = 0; i < num_rows; i++, m_data += stride)
      blas_axpy(dim, alpha, m_data, 1, data, 1);

  } else {
    Vector<Real> ones(M.NumRows());
    ones.Set(1.0);
    this->AddMatVec(alpha, M, kTrans, ones, beta);
  }
}

template<typename Real>
void VectorBase<Real>::AddDiagMat2(
    Real alpha, const MatrixBase<Real> &M,
    MatrixTransposeType trans, Real beta) {
  if (trans == kNoTrans) {
    KALDI_ASSERT(this->dim_ == M.NumRows());
    MatrixIndexT rows = this->dim_, cols = M.NumCols(),
           mat_stride = M.Stride();
    Real *data = this->data_;
    const Real *mat_data = M.Data();
    for (MatrixIndexT i = 0; i < rows; i++, mat_data += mat_stride, data++)
      *data = beta * *data + alpha * blas_dot(cols,mat_data,1,mat_data,1);
  } else {
    KALDI_ASSERT(this->dim_ == M.NumCols());
    MatrixIndexT rows = M.NumRows(), cols = this->dim_,
           mat_stride = M.Stride();
    Real *data = this->data_;
    const Real *mat_data = M.Data();
    for (MatrixIndexT i = 0; i < cols; i++, mat_data++, data++)
      *data = beta * *data + alpha * blas_dot(rows, mat_data, mat_stride,
                                                 mat_data, mat_stride);
  }
}

template<typename Real>
void VectorBase<Real>::Read(std::istream & is,  bool binary, bool add) {
  if (add) {
    Vector<Real> tmp(Dim());
    tmp.Read(is, binary, false);  // read without adding.
    if (this->Dim() != tmp.Dim()) {
      KALDI_ERR << "VectorBase::Read, size mismatch " << this->Dim()<<" vs. "<<tmp.Dim();
    }
    this->AddVec(1.0, tmp);
    return;
  } // now assume add == false.

  //  In order to avoid rewriting this, we just declare a Vector and
  // use it to read the data, then copy.
  Vector<Real> tmp;
  tmp.Read(is, binary, false);
  if (tmp.Dim() != Dim())
    KALDI_ERR << "VectorBase<Real>::Read, size mismatch "
              << Dim() << " vs. " << tmp.Dim();
  CopyFromVec(tmp);
}

template<typename Real>
void VectorBase<Real>::Write(std::ostream & os, bool binary) const {
  if (!os.good()) {
    KALDI_ERR << "Failed to write vector to stream: stream not good";
  }
  if (binary) {
    std::string my_token = (sizeof(Real) == 4 ? "FV" : "DV");
    WriteToken(os, binary, my_token);

    int32 size = Dim();  // make the size 32-bit on disk.
    KALDI_ASSERT(Dim() == (MatrixIndexT) size);
    WriteBasicType(os, binary, size);
    os.write(reinterpret_cast<const char*>(Data()), sizeof(Real) * size);
  } else {
    os << " [ ";
    for (MatrixIndexT i = 0; i < Dim(); i++)
      os << (*this)(i) << " ";
    os << "]\n";
  }
  if (!os.good())
    KALDI_ERR << "Failed to write vector to stream";
}

template<typename Real>
void Vector<Real>::Swap(Vector<Real> *other) {
  std::swap(this->data_, other->data_);
  std::swap(this->dim_, other->dim_);
}

template<typename Real>
void Vector<Real>::Read(std::istream & is,  bool binary, bool add) {
  if (add) {
    Vector<Real> tmp(this->Dim());
    tmp.Read(is, binary, false);  // read without adding.
    if (this->Dim() == 0) this->Resize(tmp.Dim());
    if (this->Dim() != tmp.Dim()) {
      KALDI_ERR << "Vector<Real>::Read, adding but dimensions mismatch "
                << this->Dim() << " vs. " << tmp.Dim();
    }
    this->AddVec(1.0, tmp);
    return;
  } // now assume add == false.

  std::ostringstream specific_error;
  MatrixIndexT pos_at_start = is.tellg();

  if (binary) {
    int peekval = Peek(is, binary);
    const char *my_token =  (sizeof(Real) == 4 ? "FV" : "DV");
    char other_token_start = (sizeof(Real) == 4 ? 'D' : 'F');
    if (peekval == other_token_start) {  // need to instantiate the other type to read it.
      typedef typename OtherReal<Real>::Real OtherType;  // if Real == float, OtherType == double, and vice versa.
      Vector<OtherType> other(this->Dim());
      other.Read(is, binary, false);  // add is false at this point.
      if (this->Dim() != other.Dim()) this->Resize(other.Dim());
      this->CopyFromVec(other);
      return;
    }
    std::string token;
    ReadToken(is, binary, &token);
    if (token != my_token) {
      specific_error << ": Expected token " << my_token << ", got " << token;
      goto bad;
    }
    int32 size;
    ReadBasicType(is, binary, &size);  // throws on error.
    if ((MatrixIndexT)size != this->Dim())  this->Resize(size);
    if (size > 0)
      is.read(reinterpret_cast<char*>(this->data_), sizeof(Real)*size);
    if (is.fail()) {
      specific_error << "Error reading vector data (binary mode); truncated "
          "stream? (size = " << size << ")";
      goto bad;
    }
    return;
  } else {  // Text mode reading; format is " [ 1.1 2.0 3.4 ]\n"
    std::string s;
    is >> s;
    // if ((s.compare("DV") == 0) || (s.compare("FV") == 0)) {  // Back compatibility.
    //  is >> s;  // get dimension
    //  is >> s;  // get "["
    // }
    if (is.fail()) { specific_error << "EOF while trying to read vector."; goto bad; }
    if (s.compare("[]") == 0) { Resize(0); return; } // tolerate this variant.
    if (s.compare("[")) { specific_error << "Expected \"[\" but got " << s; goto bad; }
    std::vector<Real> data;
    while (1) {
      int i = is.peek();
      if (i == '-' || (i >= '0' && i <= '9')) {  // common cases first.
        Real r;
        is >> r;
        if (is.fail()) { specific_error << "Failed to read number."; goto bad; }
        if (! std::isspace(is.peek()) && is.peek() != ']') {
          specific_error << "Expected whitespace after number."; goto bad;
        }
        data.push_back(r);
        // But don't eat whitespace... we want to check that it's not newlines
        // which would be valid only for a matrix.
      } else if (i == ' ' || i == '\t') {
        is.get();
      } else if (i == ']') {
        is.get();  // eat the ']'
        this->Resize(data.size());
        for (size_t j = 0; j < data.size(); j++)
          this->data_[j] = data[j];
        i = is.peek();
        if (static_cast<char>(i) == '\r') {
          is.get();
          is.get();  // get \r\n (must eat what we wrote)
        } else if (static_cast<char>(i) == '\n') { is.get(); } // get \n (must eat what we wrote)
        if (is.fail()) {
          KALDI_WARN << "After end of vector data, read error.";
          // we got the data we needed, so just warn for this error.
        }
        return;  // success.
      } else if (i == -1) {
        specific_error << "EOF while reading vector data.";
        goto bad;
      } else if (i == '\n' || i == '\r') {
        specific_error << "Newline found while reading vector (maybe it's a matrix?)";
        goto bad;
      } else {
        is >> s;  // read string.
        if (!KALDI_STRCASECMP(s.c_str(), "inf") ||
            !KALDI_STRCASECMP(s.c_str(), "infinity")) {
          data.push_back(std::numeric_limits<Real>::infinity());
          KALDI_WARN << "Reading infinite value into vector.";
        } else if (!KALDI_STRCASECMP(s.c_str(), "nan")) {
          data.push_back(std::numeric_limits<Real>::quiet_NaN());
          KALDI_WARN << "Reading NaN value into vector.";
        } else {
          specific_error << "Expecting numeric vector data, got " << s;
          goto  bad;
        }
      }
    }
  }
  // we never reach this line (the while loop returns directly).
bad:
  KALDI_ERR << "Failed to read vector from stream.  " << specific_error.str()
            << " File position at start is "
            << pos_at_start<<", currently "<<is.tellg();
}

template<typename Real>
void Vector<Real>::Resize(const MatrixIndexT dim, MatrixResizeType resize_type) {

  // the next block uses recursion to handle what we have to do if
  // resize_type == kCopyData.
  if (resize_type == kCopyData) {
    if (this->data_ == NULL || dim == 0) resize_type = kSetZero;  // nothing to copy.
    else if (this->dim_ == dim) { return; } // nothing to do.
    else {
      // set tmp to a vector of the desired size.
      Vector<Real> tmp(dim, kUndefined);
      if (dim > this->dim_) {
        memcpy(tmp.data_, this->data_, sizeof(Real)*this->dim_);
        memset(tmp.data_+this->dim_, 0, sizeof(Real)*(dim-this->dim_));
      } else {
        memcpy(tmp.data_, this->data_, sizeof(Real)*dim);
      }
      tmp.Swap(this);
      // and now let tmp go out of scope, deleting what was in *this.
      return;
    }
  }
  // At this point, resize_type == kSetZero or kUndefined.

  if (this->data_ != NULL) {
    if (this->dim_ == dim) {
      if (resize_type == kSetZero) this->SetZero();
      return;
    } else {
      Destroy();
    }
  }
  Init(dim);
  if (resize_type == kSetZero) this->SetZero();
}

template<typename Real>
inline void Vector<Real>::Init(const MatrixIndexT dim) {
  KALDI_ASSERT(dim >= 0);
  if (dim == 0) {
    this->dim_ = 0;
    this->data_ = NULL;
    return;
  }
  MatrixIndexT size;
  void *data;
  void *free_data;

  size = dim * sizeof(Real);

  if ((data = KALDI_MEMALIGN(16, size, &free_data)) != NULL) {
    this->data_ = static_cast<Real*> (data);
    this->dim_ = dim;
  } else {
    throw std::bad_alloc();
  }
}

/// Deallocates memory and sets object to empty vector.
template<typename Real>
void Vector<Real>::Destroy() {
  /// we need to free the data block if it was defined
  if (this->data_ != NULL)
    KALDI_MEMALIGN_FREE(this->data_);
  this->data_ = NULL;
  this->dim_ = 0;
}

template<typename Real>
Real VecVec(const VectorBase<Real> &a,
            const VectorBase<Real> &b) {
  MatrixIndexT adim = a.Dim();
  KALDI_ASSERT(adim == b.Dim());
  return blas_dot(adim, a.Data(), 1, b.Data(), 1);
}

template
float VecVec<>(const VectorBase<float> &a,
               const VectorBase<float> &b);
template
double VecVec<>(const VectorBase<double> &a,
                const VectorBase<double> &b);

template<typename Real, typename OtherReal>
Real VecVec(const VectorBase<Real> &ra,
            const VectorBase<OtherReal> &rb) {
  MatrixIndexT adim = ra.Dim();
  KALDI_ASSERT(adim == rb.Dim());
  const Real *a_data = ra.Data();
  const OtherReal *b_data = rb.Data();
  Real sum = 0.0;
  for (MatrixIndexT i = 0; i < adim; i++)
    sum += a_data[i]*b_data[i];
  return sum;
}

// instantiate the template above.
template
float VecVec<>(const VectorBase<float> &ra,
               const VectorBase<double> &rb);
template
double VecVec<>(const VectorBase<double> &ra,
                const VectorBase<float> &rb);


template<>
template<>
void VectorBase<float>::AddVec(const float alpha,
                               const VectorBase<float> &v) {
  KALDI_ASSERT(dim_ == v.dim_);
  KALDI_ASSERT(&v != this);
  blas_axpy(dim_, alpha, v.Data(), 1, data_, 1);
}

template<>
template<>
void VectorBase<double>::AddVec(const double alpha,
                                const VectorBase<double> &v) {
  KALDI_ASSERT(dim_ == v.dim_);
  KALDI_ASSERT(&v != this);
  blas_axpy(dim_, alpha, v.Data(), 1, data_, 1);
}

template class Vector<float>;
template class Vector<double>;
template class VectorBase<float>;
template class VectorBase<double>;

}  // namespace kaldi
