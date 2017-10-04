// matrix/kaldi-matrix.cc

// Copyright 2009-2011   Lukas Burget;  Ondrej Glembek;  Go Vivace Inc.;
//                       Microsoft Corporation;  Saarland University;
//                       Yanmin Qian;  Petr Schwarz;  Jan Silovsky;
//                       Haihua Xu

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

#include "matrix/kaldi-matrix.h"
#include "matrix/cblas-wrappers.h"

namespace kaldi {

template<typename Real>
void MatrixBase<Real>::SetZero() {
  if (num_cols_ == stride_)
    memset(data_, 0, sizeof(Real)*num_rows_*num_cols_);
  else
    for (MatrixIndexT row = 0; row < num_rows_; row++)
      memset(data_ + row*stride_, 0, sizeof(Real)*num_cols_);
}

template<typename Real>
void MatrixBase<Real>::Set(Real value) {
  for (MatrixIndexT row = 0; row < num_rows_; row++) {
    for (MatrixIndexT col = 0; col < num_cols_; col++) {
      (*this)(row, col) = value;
    }
  }
}

template<typename Real>
void MatrixBase<Real>::SetRandn() {
  kaldi::RandomState rstate;
  for (MatrixIndexT row = 0; row < num_rows_; row++) {
    Real *row_data = this->RowData(row);
    MatrixIndexT nc = (num_cols_ % 2 == 1) ? num_cols_ - 1 : num_cols_;
    for (MatrixIndexT col = 0; col < nc; col += 2) {
      kaldi::RandGauss2(row_data + col, row_data + col + 1, &rstate);
    }
    if (nc != num_cols_) row_data[nc] = static_cast<Real>(kaldi::RandGauss(&rstate));
  }
}

template<typename Real>
template<typename OtherReal>
void MatrixBase<Real>::CopyFromMat(const MatrixBase<OtherReal> &M,
                                   MatrixTransposeType Trans) {
  if (sizeof(Real) == sizeof(OtherReal) &&
      static_cast<const void*>(M.Data()) ==
      static_cast<const void*>(this->Data())) {
    // CopyFromMat called on same data.  Nothing to do (except sanity checks).
    KALDI_ASSERT(Trans == kNoTrans && M.NumRows() == NumRows() &&
                 M.NumCols() == NumCols() && M.Stride() == Stride());
    return;
  }
  if (Trans == kNoTrans) {
    KALDI_ASSERT(num_rows_ == M.NumRows() && num_cols_ == M.NumCols());
    for (MatrixIndexT i = 0; i < num_rows_; i++)
      (*this).Row(i).CopyFromVec(M.Row(i));
  } else {
    KALDI_ASSERT(num_cols_ == M.NumRows() && num_rows_ == M.NumCols());
    int32 this_stride = stride_, other_stride = M.Stride();
    Real *this_data = data_;
    const OtherReal *other_data = M.Data();
    for (MatrixIndexT i = 0; i < num_rows_; i++)
      for (MatrixIndexT j = 0; j < num_cols_; j++)
        this_data[i * this_stride + j] = other_data[j * other_stride + i];
  }
}

// template instantiations.
template
void MatrixBase<float>::CopyFromMat(const MatrixBase<double> & M,
                                    MatrixTransposeType Trans);
template
void MatrixBase<double>::CopyFromMat(const MatrixBase<float> & M,
                                     MatrixTransposeType Trans);
template
void MatrixBase<float>::CopyFromMat(const MatrixBase<float> & M,
                                    MatrixTransposeType Trans);
template
void MatrixBase<double>::CopyFromMat(const MatrixBase<double> & M,
                                     MatrixTransposeType Trans);
 
template<typename Real>
void MatrixBase<Real>::CopyColFromVec(const VectorBase<Real> &rv,
                                      const MatrixIndexT col) {
  KALDI_ASSERT(rv.Dim() == num_rows_ &&
               static_cast<UnsignedMatrixIndexT>(col) <
               static_cast<UnsignedMatrixIndexT>(num_cols_));

  const Real *rv_data = rv.Data();
  Real *col_data = data_ + col;

  for (MatrixIndexT r = 0; r < num_rows_; r++)
    col_data[r * stride_] = rv_data[r];
}

template<typename Real>
Real MatrixBase<Real>::Trace(bool check_square) const  {
  KALDI_ASSERT(!check_square || num_rows_ == num_cols_);
  Real ans = 0.0;
  for (MatrixIndexT r = 0;r < std::min(num_rows_, num_cols_);r++) ans += data_ [r + stride_*r];
  return ans;
}

template<typename Real>
Real MatrixBase<Real>::Max() const {
  KALDI_ASSERT(num_rows_ > 0 && num_cols_ > 0);
  Real ans= *data_;
  for (MatrixIndexT r = 0; r < num_rows_; r++)
    for (MatrixIndexT c = 0; c < num_cols_; c++)
      if (data_[c + stride_*r] > ans)
        ans = data_[c + stride_*r];
  return ans;
}

template<typename Real>
Real MatrixBase<Real>::Min() const {
  KALDI_ASSERT(num_rows_ > 0 && num_cols_ > 0);
  Real ans= *data_;
  for (MatrixIndexT r = 0; r < num_rows_; r++)
    for (MatrixIndexT c = 0; c < num_cols_; c++)
      if (data_[c + stride_*r] < ans)
        ans = data_[c + stride_*r];
  return ans;
}

template<typename Real>
void MatrixBase<Real>::MulElements(const MatrixBase<Real> &a) {
  KALDI_ASSERT(a.NumRows() == num_rows_ && a.NumCols() == num_cols_);

  if (num_cols_ == stride_ && num_cols_ == a.stride_) {
    mul_elements(num_rows_ * num_cols_, a.data_, data_);
  } else {
    MatrixIndexT a_stride = a.stride_, stride = stride_;
    Real *data = data_, *a_data = a.data_;
    for (MatrixIndexT i = 0; i < num_rows_; i++) {
      mul_elements(num_cols_, a_data, data);
      a_data += a_stride;
      data += stride;
    }
  }
}

template<typename Real> void MatrixBase<Real>::Scale(Real alpha) {
  if (alpha == 1.0) return;
  if (num_rows_ == 0) return;
  blas_scale(data_, alpha, num_rows_, num_cols_, stride_);
}

template<typename Real> void MatrixBase<Real>::Max(const MatrixBase<Real> &A) {
  KALDI_ASSERT(A.NumRows() == NumRows() && A.NumCols() == NumCols());
  for (MatrixIndexT row = 0; row < num_rows_; row++) {
    Real *row_data = RowData(row);
    const Real *other_row_data = A.RowData(row);
    MatrixIndexT num_cols = num_cols_;
    for (MatrixIndexT col = 0; col < num_cols; col++) {
      row_data[col] = std::max(row_data[col],
                               other_row_data[col]);
    }
  }
}

template<typename Real> void MatrixBase<Real>::Min(const MatrixBase<Real> &A) {
  KALDI_ASSERT(A.NumRows() == NumRows() && A.NumCols() == NumCols());
  for (MatrixIndexT row = 0; row < num_rows_; row++) {
    Real *row_data = RowData(row);
    const Real *other_row_data = A.RowData(row);
    MatrixIndexT num_cols = num_cols_;
    for (MatrixIndexT col = 0; col < num_cols; col++) {
      row_data[col] = std::min(row_data[col],
                               other_row_data[col]);
    }
  }
}

template<typename Real>
void MatrixBase<Real>::Transpose() {
  KALDI_ASSERT(num_rows_ == num_cols_);
  MatrixIndexT M = num_rows_;
  for (MatrixIndexT i = 0;i < M;i++) {
    for (MatrixIndexT j = 0;j < i;j++) {
      Real &a = (*this)(i, j), &b = (*this)(j, i);
      std::swap(a, b);
    }
  }
}

template<typename Real>
void MatrixBase<Real>::ApplyFloor(Real floor_val) {
  MatrixIndexT num_rows = num_rows_, num_cols = num_cols_;
  for (MatrixIndexT i = 0; i < num_rows; i++) {
    Real *data = this->RowData(i);
    for (MatrixIndexT j = 0; j < num_cols; j++)
      data[j] = (data[j] < floor_val ? floor_val : data[j]);
  }
}

template<typename Real>
bool MatrixBase<Real>::IsZero(Real cutoff)const {
  MatrixIndexT R = num_rows_, C = num_cols_;
  Real bad_max = 0.0;
  for (MatrixIndexT i = 0;i < R;i++)
    for (MatrixIndexT j = 0;j < C;j++)
      bad_max = std::max(bad_max, static_cast<Real>(std::abs( (*this)(i, j) )));
  return (bad_max <= cutoff);
}

template<typename Real>
bool MatrixBase<Real>::Equal(const MatrixBase<Real> &other) const {
  if (num_rows_ != other.num_rows_ || num_cols_ != other.num_cols_)
    KALDI_ERR << "Equal: size mismatch.";
  for (MatrixIndexT i = 0; i < num_rows_; i++)
    for (MatrixIndexT j = 0; j < num_cols_; j++)
      if ( (*this)(i, j) != other(i, j))
        return false;
  return true;
}

template<typename Real>
void MatrixBase<Real>::Add(const Real alpha) {
  Real *data = data_;
  MatrixIndexT stride = stride_;
  for (MatrixIndexT r = 0; r < num_rows_; r++)
    for (MatrixIndexT c = 0; c < num_cols_; c++)
      data[c + stride*r] += alpha;
}

template<typename Real>
void MatrixBase<Real>::AddMat(const Real alpha, const MatrixBase<Real>& A,
                               MatrixTransposeType transA) {
  if (&A == this) {
    if (transA == kNoTrans) {
      Scale(alpha + 1.0);
    } else {
      KALDI_ASSERT(num_rows_ == num_cols_ && "AddMat: adding to self (transposed): not symmetric.");
      Real *data = data_;
      if (alpha == 1.0) {  // common case-- handle separately.
        for (MatrixIndexT row = 0; row < num_rows_; row++) {
          for (MatrixIndexT col = 0; col < row; col++) {
            Real *lower = data + (row * stride_) + col, *upper = data + (col
                                                                          * stride_) + row;
            Real sum = *lower + *upper;
            *lower = *upper = sum;
          }
          *(data + (row * stride_) + row) *= 2.0;  // diagonal.
        }
      } else {
        for (MatrixIndexT row = 0; row < num_rows_; row++) {
          for (MatrixIndexT col = 0; col < row; col++) {
            Real *lower = data + (row * stride_) + col, *upper = data + (col
                                                                          * stride_) + row;
            Real lower_tmp = *lower;
            *lower += alpha * *upper;
            *upper += alpha * lower_tmp;
          }
          *(data + (row * stride_) + row) *= (1.0 + alpha);  // diagonal.
        }
      }
    }
  } else {
    int aStride = (int) A.stride_, stride = stride_;
    Real *adata = A.data_, *data = data_;
    if (transA == kNoTrans) {
      KALDI_ASSERT(A.num_rows_ == num_rows_ && A.num_cols_ == num_cols_);
      if (num_rows_ == 0) return;
      for (MatrixIndexT row = 0; row < num_rows_; row++, adata += aStride,
               data += stride) {
        blas_axpy(num_cols_, alpha, adata, 1, data, 1);
      }
    } else {
      KALDI_ASSERT(A.num_cols_ == num_rows_ && A.num_rows_ == num_cols_);
      if (num_rows_ == 0) return;
      for (MatrixIndexT row = 0; row < num_rows_; row++, adata++, data += stride)
        blas_axpy(num_cols_, alpha, adata, aStride, data, 1);
    }
  }
}

template<typename Real>
void MatrixBase<Real>::Read(std::istream & is, bool binary, bool add) {
  if (add) {
    Matrix<Real> tmp(num_rows_, num_cols_);
    tmp.Read(is, binary, false);  // read without adding.
    if (tmp.num_rows_ != this->num_rows_ || tmp.num_cols_ != this->num_cols_)
      KALDI_ERR << "MatrixBase::Read, size mismatch "
                << this->num_rows_ << ", " << this->num_cols_
                << " vs. " << tmp.num_rows_ << ", " << tmp.num_cols_;
    this->AddMat(1.0, tmp);
    return;
  }
  // now assume add == false.

  //  In order to avoid rewriting this, we just declare a Matrix and
  // use it to read the data, then copy.
  Matrix<Real> tmp;
  tmp.Read(is, binary, false);
  if (tmp.NumRows() != NumRows() || tmp.NumCols() != NumCols()) {
    KALDI_ERR << "MatrixBase<Real>::Read, size mismatch "
              << NumRows() << " x " << NumCols() << " versus "
              << tmp.NumRows() << " x " << tmp.NumCols();
  }
  CopyFromMat(tmp);
}

template<typename Real>
void MatrixBase<Real>::Write(std::ostream &os, bool binary) const {
  if (!os.good()) {
    KALDI_ERR << "Failed to write matrix to stream: stream not good";
  }
  if (binary) {  // Use separate binary and text formats,
    // since in binary mode we need to know if it's float or double.
    std::string my_token = (sizeof(Real) == 4 ? "FM" : "DM");

    WriteToken(os, binary, my_token);
    {
      int32 rows = this->num_rows_;  // make the size 32-bit on disk.
      int32 cols = this->num_cols_;
      KALDI_ASSERT(this->num_rows_ == (MatrixIndexT) rows);
      KALDI_ASSERT(this->num_cols_ == (MatrixIndexT) cols);
      WriteBasicType(os, binary, rows);
      WriteBasicType(os, binary, cols);
    }
    if (Stride() == NumCols())
      os.write(reinterpret_cast<const char*> (Data()), sizeof(Real)
               * static_cast<size_t>(num_rows_) * static_cast<size_t>(num_cols_));
    else
      for (MatrixIndexT i = 0; i < num_rows_; i++)
        os.write(reinterpret_cast<const char*> (RowData(i)), sizeof(Real)
                 * num_cols_);
    if (!os.good()) {
      KALDI_ERR << "Failed to write matrix to stream";
    }
  } else {  // text mode.
    if (num_cols_ == 0) {
      os << " [ ]\n";
    } else {
      os << " [";
      for (MatrixIndexT i = 0; i < num_rows_; i++) {
        os << "\n  ";
        for (MatrixIndexT j = 0; j < num_cols_; j++)
          os << (*this)(i, j) << " ";
      }
      os << "]\n";
    }
  }
}

template<typename Real>
void Matrix<Real>::Swap(Matrix<Real> *other) {
  std::swap(this->data_, other->data_);
  std::swap(this->num_cols_, other->num_cols_);
  std::swap(this->num_rows_, other->num_rows_);
  std::swap(this->stride_, other->stride_);
}

// Copy constructor.  Copies data to newly allocated memory.
template<typename Real>
Matrix<Real>::Matrix (const MatrixBase<Real> & M,
                      MatrixTransposeType trans/*=kNoTrans*/)
    : MatrixBase<Real>() {
  if (trans == kNoTrans) {
    Resize(M.num_rows_, M.num_cols_);
    this->CopyFromMat(M);
  } else {
    Resize(M.num_cols_, M.num_rows_);
    this->CopyFromMat(M, kTrans);
  }
}

// Copy constructor.  Copies data to newly allocated memory.
template<typename Real>
Matrix<Real>::Matrix (const Matrix<Real> & M):
    MatrixBase<Real>() {
  Resize(M.num_rows_, M.num_cols_);
  this->CopyFromMat(M);
}

/// Copy constructor from another type.
template<typename Real>
template<typename OtherReal>
Matrix<Real>::Matrix(const MatrixBase<OtherReal> & M,
                     MatrixTransposeType trans) : MatrixBase<Real>() {
  if (trans == kNoTrans) {
    Resize(M.NumRows(), M.NumCols());
    this->CopyFromMat(M);
  } else {
    Resize(M.NumCols(), M.NumRows());
    this->CopyFromMat(M, kTrans);
  }
}

// Instantiate this constructor for float->double and double->float.
template
Matrix<float>::Matrix(const MatrixBase<double> & M,
                      MatrixTransposeType trans);
template
Matrix<double>::Matrix(const MatrixBase<float> & M,
                       MatrixTransposeType trans);

template<typename Real>
void Matrix<Real>::Read(std::istream & is, bool binary, bool add) {
  if (add) {
    Matrix<Real> tmp;
    tmp.Read(is, binary, false);  // read without adding.
    if (this->num_rows_ == 0) this->Resize(tmp.num_rows_, tmp.num_cols_);
    else {
      if (this->num_rows_ != tmp.num_rows_ || this->num_cols_ != tmp.num_cols_) {
        if (tmp.num_rows_ == 0) return;  // do nothing in this case.
        else KALDI_ERR << "Matrix::Read, size mismatch "
                       << this->num_rows_ <<  ", " << this->num_cols_
                       << " vs. " << tmp.num_rows_ << ", " << tmp.num_cols_;
      }
    }
    this->AddMat(1.0, tmp);
    return;
  }

  // now assume add == false.
  MatrixIndexT pos_at_start = is.tellg();
  std::ostringstream specific_error;

  if (binary) {  // Read in binary mode.
    int peekval = Peek(is, binary);
    const char *my_token =  (sizeof(Real) == 4 ? "FM" : "DM");
    char other_token_start = (sizeof(Real) == 4 ? 'D' : 'F');
    if (peekval == other_token_start) {  // need to instantiate the other type to read it.
      typedef typename OtherReal<Real>::Real OtherType;  // if Real == float, OtherType == double, and vice versa.
      Matrix<OtherType> other(this->num_rows_, this->num_cols_);
      other.Read(is, binary, false);  // add is false at this point anyway.
      this->Resize(other.NumRows(), other.NumCols());
      this->CopyFromMat(other);
      return;
    }
    std::string token;
    ReadToken(is, binary, &token);
    if (token != my_token) {
      specific_error << ": Expected token " << my_token << ", got " << token;
      goto bad;
    }
    int32 rows, cols;
    ReadBasicType(is, binary, &rows);  // throws on error.
    ReadBasicType(is, binary, &cols);  // throws on error.
    if ((MatrixIndexT)rows != this->num_rows_ || (MatrixIndexT)cols != this->num_cols_) {
      this->Resize(rows, cols);
    }
    if (this->Stride() == this->NumCols() && rows*cols!=0) {
      is.read(reinterpret_cast<char*>(this->Data()),
              sizeof(Real)*rows*cols);
      if (is.fail()) goto bad;
    } else {
      for (MatrixIndexT i = 0; i < (MatrixIndexT)rows; i++) {
        is.read(reinterpret_cast<char*>(this->RowData(i)), sizeof(Real)*cols);
        if (is.fail()) goto bad;
      }
    }
    if (is.eof()) return;
    if (is.fail()) goto bad;
    return;
  } else {  // Text mode.
    std::string str;
    is >> str; // get a token
    if (is.fail()) { specific_error << ": Expected \"[\", got EOF"; goto bad; }
    // if ((str.compare("DM") == 0) || (str.compare("FM") == 0)) {  // Back compatibility.
    // is >> str;  // get #rows
    //  is >> str;  // get #cols
    //  is >> str;  // get "["
    // }
    if (str == "[]") { Resize(0, 0); return; } // Be tolerant of variants.
    else if (str != "[") {
      specific_error << ": Expected \"[\", got \"" << str << '"';
      goto bad;
    }
    // At this point, we have read "[".
    std::vector<std::vector<Real>* > data;
    std::vector<Real> *cur_row = new std::vector<Real>;
    while (1) {
      int i = is.peek();
      if (i == -1) { specific_error << "Got EOF while reading matrix data"; goto cleanup; }
      else if (static_cast<char>(i) == ']') {  // Finished reading matrix.
        is.get();  // eat the "]".
        i = is.peek();
        if (static_cast<char>(i) == '\r') {
          is.get();
          is.get();  // get \r\n (must eat what we wrote)
        } else if (static_cast<char>(i) == '\n') { is.get(); } // get \n (must eat what we wrote)
        if (is.fail()) {
          KALDI_WARN << "After end of matrix data, read error.";
          // we got the data we needed, so just warn for this error.
        }
        // Now process the data.
        if (!cur_row->empty()) data.push_back(cur_row);
        else delete(cur_row);
        cur_row = NULL;
        if (data.empty()) { this->Resize(0, 0); return; }
        else {
          int32 num_rows = data.size(), num_cols = data[0]->size();
          this->Resize(num_rows, num_cols);
          for (int32 i = 0; i < num_rows; i++) {
            if (static_cast<int32>(data[i]->size()) != num_cols) {
              specific_error << "Matrix has inconsistent #cols: " << num_cols
                             << " vs." << data[i]->size() << " (processing row"
                             << i << ")";
              goto cleanup;
            }
            for (int32 j = 0; j < num_cols; j++)
              (*this)(i, j) = (*(data[i]))[j];
            delete data[i];
            data[i] = NULL;
          }
        }
        return;
      } else if (static_cast<char>(i) == '\n' || static_cast<char>(i) == ';') {
        // End of matrix row.
        is.get();
        if (cur_row->size() != 0) {
          data.push_back(cur_row);
          cur_row = new std::vector<Real>;
          cur_row->reserve(data.back()->size());
        }
      } else if ( (i >= '0' && i <= '9') || i == '-' ) {  // A number...
        Real r;
        is >> r;
        if (is.fail()) {
          specific_error << "Stream failure/EOF while reading matrix data.";
          goto cleanup;
        }
        cur_row->push_back(r);
      } else if (isspace(i)) {
        is.get();  // eat the space and do nothing.
      } else {  // NaN or inf or error.
        std::string str;
        is >> str;
        if (!KALDI_STRCASECMP(str.c_str(), "inf") ||
            !KALDI_STRCASECMP(str.c_str(), "infinity")) {
          cur_row->push_back(std::numeric_limits<Real>::infinity());
          KALDI_WARN << "Reading infinite value into matrix.";
        } else if (!KALDI_STRCASECMP(str.c_str(), "nan")) {
          cur_row->push_back(std::numeric_limits<Real>::quiet_NaN());
          KALDI_WARN << "Reading NaN value into matrix.";
        } else {
          specific_error << "Expecting numeric matrix data, got " << str;
          goto cleanup;
        }
      }
    }
    // Note, we never leave the while () loop before this
    // line (we return from it.)
 cleanup: // We only reach here in case of error in the while loop above.
    if(cur_row != NULL)
      delete cur_row;
    for (size_t i = 0; i < data.size(); i++)
      if(data[i] != NULL)
        delete data[i];
    // and then go on to "bad" below, where we print error.
  }
bad:
  KALDI_ERR << "Failed to read matrix from stream.  " << specific_error.str()
            << " File position at start is "
            << pos_at_start << ", currently " << is.tellg();
}

template<typename Real>
void Matrix<Real>::Transpose() {
  if (this->num_rows_ != this->num_cols_) {
    Matrix<Real> tmp(*this, kTrans);
    Resize(this->num_cols_, this->num_rows_);
    this->CopyFromMat(tmp);
  } else {
    (static_cast<MatrixBase<Real>&>(*this)).Transpose();
  }
}

template<typename Real>
void Matrix<Real>::Resize(const MatrixIndexT rows,
                          const MatrixIndexT cols,
                          MatrixResizeType resize_type,
                          MatrixStrideType stride_type) {
  // the next block uses recursion to handle what we have to do if
  // resize_type == kCopyData.
  if (resize_type == kCopyData) {
    if (this->data_ == NULL || rows == 0) resize_type = kSetZero;  // nothing to copy.
    else if (rows == this->num_rows_ && cols == this->num_cols_) { return; } // nothing to do.
    else {
      // set tmp to a matrix of the desired size; if new matrix
      // is bigger in some dimension, zero it.
      MatrixResizeType new_resize_type =
          (rows > this->num_rows_ || cols > this->num_cols_) ? kSetZero : kUndefined;
      Matrix<Real> tmp(rows, cols, new_resize_type);
      MatrixIndexT rows_min = std::min(rows, this->num_rows_),
          cols_min = std::min(cols, this->num_cols_);
      tmp.Range(0, rows_min, 0, cols_min).
          CopyFromMat(this->Range(0, rows_min, 0, cols_min));
      tmp.Swap(this);
      // and now let tmp go out of scope, deleting what was in *this.
      return;
    }
  }
  // At this point, resize_type == kSetZero or kUndefined.

  if (MatrixBase<Real>::data_ != NULL) {
    if (rows == MatrixBase<Real>::num_rows_
        && cols == MatrixBase<Real>::num_cols_) {
      if (resize_type == kSetZero)
        this->SetZero();
      return;
    }
    else
      Destroy();
  }
  Init(rows, cols, stride_type);
  if (resize_type == kSetZero) MatrixBase<Real>::SetZero();
}

template<typename Real>
void Matrix<Real>::Destroy() {
  // we need to free the data block if it was defined
  if (NULL != MatrixBase<Real>::data_)
    KALDI_MEMALIGN_FREE( MatrixBase<Real>::data_);
  MatrixBase<Real>::data_ = NULL;
  MatrixBase<Real>::num_rows_ = MatrixBase<Real>::num_cols_
      = MatrixBase<Real>::stride_ = 0;
}

template<typename Real>
inline void Matrix<Real>::Init(const MatrixIndexT rows,
                               const MatrixIndexT cols,
                               const MatrixStrideType stride_type) {
  if (rows * cols == 0) {
    KALDI_ASSERT(rows == 0 && cols == 0);
    this->num_rows_ = 0;
    this->num_cols_ = 0;
    this->stride_ = 0;
    this->data_ = NULL;
    return;
  }
  KALDI_ASSERT(rows > 0 && cols > 0);
  MatrixIndexT skip, stride;
  size_t size;
  void *data;  // aligned memory block
  void *temp;  // memory block to be really freed

  // compute the size of skip and real cols
  skip = ((16 / sizeof(Real)) - cols % (16 / sizeof(Real)))
      % (16 / sizeof(Real));
  stride = cols + skip;
  size = static_cast<size_t>(rows) * static_cast<size_t>(stride)
      * sizeof(Real);

  // allocate the memory and set the right dimensions and parameters
  if (NULL != (data = KALDI_MEMALIGN(16, size, &temp))) {
    MatrixBase<Real>::data_        = static_cast<Real *> (data);
    MatrixBase<Real>::num_rows_      = rows;
    MatrixBase<Real>::num_cols_      = cols;
    MatrixBase<Real>::stride_  = (stride_type == kDefaultStride ? stride : cols);
  } else {
    throw std::bad_alloc();
  }
}

// Constructor... note that this is not const-safe as it would
// be quite complicated to implement a "const SubMatrix" class that
// would not allow its contents to be changed.
template<typename Real>
SubMatrix<Real>::SubMatrix(const MatrixBase<Real> &M,
                           const MatrixIndexT ro,
                           const MatrixIndexT r,
                           const MatrixIndexT co,
                           const MatrixIndexT c) {
  if (r == 0 || c == 0) {
    // we support the empty sub-matrix as a special case.
    KALDI_ASSERT(c == 0 && r == 0);
    this->data_ = NULL;
    this->num_cols_ = 0;
    this->num_rows_ = 0;
    this->stride_ = 0;
    return;
  }
  KALDI_ASSERT(static_cast<UnsignedMatrixIndexT>(ro) <
               static_cast<UnsignedMatrixIndexT>(M.num_rows_) &&
               static_cast<UnsignedMatrixIndexT>(co) <
               static_cast<UnsignedMatrixIndexT>(M.num_cols_) &&
               static_cast<UnsignedMatrixIndexT>(r) <=
               static_cast<UnsignedMatrixIndexT>(M.num_rows_ - ro) &&
               static_cast<UnsignedMatrixIndexT>(c) <=
               static_cast<UnsignedMatrixIndexT>(M.num_cols_ - co));
  // point to the begining of window
  MatrixBase<Real>::num_rows_ = r;
  MatrixBase<Real>::num_cols_ = c;
  MatrixBase<Real>::stride_ = M.Stride();
  MatrixBase<Real>::data_ = M.Data_workaround() +
      static_cast<size_t>(co) +
      static_cast<size_t>(ro) * static_cast<size_t>(M.Stride());
}


template<typename Real>
SubMatrix<Real>::SubMatrix(Real *data,
                           MatrixIndexT num_rows,
                           MatrixIndexT num_cols,
                           MatrixIndexT stride):
    MatrixBase<Real>(data, num_cols, num_rows, stride) { // caution: reversed order!
  if (data == NULL) {
    KALDI_ASSERT(num_rows * num_cols == 0);
    this->num_rows_ = 0;
    this->num_cols_ = 0;
    this->stride_ = 0;
  } else {
    KALDI_ASSERT(this->stride_ >= this->num_cols_);
  }
}

//Explicit instantiation of the classes
//Apparently, it seems to be necessary that the instantiation
//happens at the end of the file. Otherwise, not all the member
//functions will get instantiated.

template class Matrix<float>;
template class Matrix<double>;
template class MatrixBase<float>;
template class MatrixBase<double>;
template class SubMatrix<float>;
template class SubMatrix<double>;

} // namespace kaldi
