// matrix/kaldi-matrix.h

// Copyright 2009-2011  Ondrej Glembek;  Microsoft Corporation;  Lukas Burget;
//                      Saarland University;  Petr Schwarz;  Yanmin Qian;
//                      Karel Vesely;  Go Vivace Inc.;  Haihua Xu

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

#ifndef KALDI_MATRIX_KALDI_MATRIX_H_
#define KALDI_MATRIX_KALDI_MATRIX_H_ 1

#include "matrix/matrix-common.h"

namespace kaldi {

/// \addtogroup matrix_group
/// @{

/// Base class which provides matrix operations not involving resizing
/// or allocation.   Classes Matrix and SubMatrix inherit from it and take care
/// of allocation and resizing.
template<typename Real>
class MatrixBase {
 public:
  // so this child can access protected members of other instances.
  friend class Matrix<Real>;
 
  /// Returns number of rows (or zero for emtpy matrix).
  inline MatrixIndexT  NumRows() const { return num_rows_; }

  /// Returns number of columns (or zero for emtpy matrix).
  inline MatrixIndexT NumCols() const { return num_cols_; }

  /// Stride (distance in memory between each row).  Will be >= NumCols.
  inline MatrixIndexT Stride() const {  return stride_; }

  /// Returns size in bytes of the data held by the matrix.
  size_t  SizeInBytes() const {
    return static_cast<size_t>(num_rows_) * static_cast<size_t>(stride_) *
        sizeof(Real);
  }

  /// Gives pointer to raw data (const).
  inline const Real* Data() const {
    return data_;
  }

  /// Gives pointer to raw data (non-const).
  inline Real* Data() { return data_; }

  /// Returns pointer to data for one row (non-const)
  inline  Real* RowData(MatrixIndexT i) {
    KALDI_ASSERT(static_cast<UnsignedMatrixIndexT>(i) <
                 static_cast<UnsignedMatrixIndexT>(num_rows_));
    return data_ + i * stride_;
  }

  /// Returns pointer to data for one row (const)
  inline const Real* RowData(MatrixIndexT i) const {
    KALDI_ASSERT(static_cast<UnsignedMatrixIndexT>(i) <
                 static_cast<UnsignedMatrixIndexT>(num_rows_));
    return data_ + i * stride_;
  }

  /// Indexing operator, non-const
  /// (only checks sizes if compiled with -DKALDI_PARANOID)
  inline Real&  operator() (MatrixIndexT r, MatrixIndexT c) {
    KALDI_PARANOID_ASSERT(static_cast<UnsignedMatrixIndexT>(r) <
                          static_cast<UnsignedMatrixIndexT>(num_rows_) &&
                          static_cast<UnsignedMatrixIndexT>(c) <
                          static_cast<UnsignedMatrixIndexT>(num_cols_));
    return *(data_ + r * stride_ + c);
  }
  /// Indexing operator, provided for ease of debugging (gdb doesn't work
  /// with parenthesis operator).
  Real &Index (MatrixIndexT r, MatrixIndexT c) {  return (*this)(r, c); }

  /// Indexing operator, const
  /// (only checks sizes if compiled with -DKALDI_PARANOID)
  inline const Real operator() (MatrixIndexT r, MatrixIndexT c) const {
    KALDI_PARANOID_ASSERT(static_cast<UnsignedMatrixIndexT>(r) <
                          static_cast<UnsignedMatrixIndexT>(num_rows_) &&
                          static_cast<UnsignedMatrixIndexT>(c) <
                          static_cast<UnsignedMatrixIndexT>(num_cols_));
    return *(data_ + r * stride_ + c);
  }

  /*   Basic setting-to-special values functions. */

  /// Sets matrix to zero.
  void SetZero();
  /// Sets all elements to a specific value.
  void Set(Real);
  /// Sets to random values of a normal distribution
  void SetRandn();

  /*  Copying functions.  These do not resize the matrix! */


  /// Copy given matrix. (no resize is done).
  template<typename OtherReal>
  void CopyFromMat(const MatrixBase<OtherReal> & M,
                   MatrixTransposeType trans = kNoTrans);

  /// Copy vector into specific column of matrix.
  void CopyColFromVec(const VectorBase<Real> &v, const MatrixIndexT col);

  /* Accessing of sub-parts of the matrix. */

  /// Return specific row of matrix [const].
  inline const SubVector<Real> Row(MatrixIndexT i) const {
    KALDI_ASSERT(static_cast<UnsignedMatrixIndexT>(i) <
                 static_cast<UnsignedMatrixIndexT>(num_rows_));
    return SubVector<Real>(data_ + (i * stride_), NumCols());
  }

  /// Return specific row of matrix.
  inline SubVector<Real> Row(MatrixIndexT i) {
    KALDI_ASSERT(static_cast<UnsignedMatrixIndexT>(i) <
                 static_cast<UnsignedMatrixIndexT>(num_rows_));
    return SubVector<Real>(data_ + (i * stride_), NumCols());
  }

  /// Return a sub-part of matrix.
  inline SubMatrix<Real> Range(const MatrixIndexT row_offset,
                               const MatrixIndexT num_rows,
                               const MatrixIndexT col_offset,
                               const MatrixIndexT num_cols) const {
    return SubMatrix<Real>(*this, row_offset, num_rows,
                           col_offset, num_cols);
  }

  inline SubMatrix<Real> RowRange(const MatrixIndexT row_offset,
                                  const MatrixIndexT num_rows) const {
    return SubMatrix<Real>(*this, row_offset, num_rows, 0, num_cols_);
  }

  inline SubMatrix<Real> ColRange(const MatrixIndexT col_offset,
                                  const MatrixIndexT num_cols) const {
    return SubMatrix<Real>(*this, 0, num_rows_, col_offset, num_cols);
  }

  /* Various special functions. */
  /// Returns sum of all elements in matrix.
  Real Sum() const;
  /// Returns trace of matrix.
  Real Trace(bool check_square = true) const;
  // If check_square = true, will crash if matrix is not square.

  /// Returns maximum element of matrix.
  Real Max() const;
  /// Returns minimum element of matrix.
  Real Min() const;

  /// Element by element multiplication with a given matrix.
  void MulElements(const MatrixBase<Real> &A);

  /// Multiply each element with a scalar value.
  void Scale(Real alpha);

  /// Set, element-by-element, *this = max(*this, A)
  void Max(const MatrixBase<Real> &A);
  /// Set, element-by-element, *this = min(*this, A)
  void Min(const MatrixBase<Real> &A);

  /// Transpose the matrix.  This one is only
  /// applicable to square matrices (the one in the
  /// Matrix child class works also for non-square.
  void Transpose();

  /// Applies floor to all matrix elements
  void ApplyFloor(Real floor_val);

  /// Returns true if matrix is all zeros.
  bool IsZero(Real cutoff = 1.0e-05) const;     // replace magic number

  /// Tests for exact equality.  It's usually preferable to use ApproxEqual.
  bool Equal(const MatrixBase<Real> &other) const;

  // so it can get around const restrictions on the pointer to data_.
  friend class SubMatrix<Real>;

  /// Add a scalar to each element
  void Add(const Real alpha);

  /// *this += alpha * M [or M^T]
  void AddMat(const Real alpha, const MatrixBase<Real> &M,
              MatrixTransposeType transA = kNoTrans);

  /// stream read.
  /// Use instead of stream<<*this, if you want to add to existing contents.
  // Will throw exception on failure.
  void Read(std::istream & in, bool binary, bool add = false);
  /// write to stream.
  void Write(std::ostream & out, bool binary) const;

 protected:

  ///  Initializer, callable only from child.
  explicit MatrixBase(Real *data, MatrixIndexT cols, MatrixIndexT rows, MatrixIndexT stride) :
    data_(data), num_cols_(cols), num_rows_(rows), stride_(stride) {
    KALDI_ASSERT_IS_FLOATING_TYPE(Real);
  }

  ///  Initializer, callable only from child.
  /// Empty initializer, for un-initialized matrix.
  explicit MatrixBase(): data_(NULL) {
    KALDI_ASSERT_IS_FLOATING_TYPE(Real);
  }

  // Make sure pointers to MatrixBase cannot be deleted.
  ~MatrixBase() { }

  /// A workaround that allows SubMatrix to get a pointer to non-const data
  /// for const Matrix. Unfortunately C++ does not allow us to declare a
  /// "public const" inheritance or anything like that, so it would require
  /// a lot of work to make the SubMatrix class totally const-correct--
  /// we would have to override many of the Matrix functions.
  inline Real*  Data_workaround() const {
    return data_;
  }

  /// data memory area
  Real*   data_;

  /// these atributes store the real matrix size as it is stored in memory
  /// including memalignment
  MatrixIndexT    num_cols_;   /// < Number of columns
  MatrixIndexT    num_rows_;   /// < Number of rows
  /** True number of columns for the internal matrix. This number may differ
   * from num_cols_ as memory alignment might be used. */
  MatrixIndexT    stride_;
 private:
  KALDI_DISALLOW_COPY_AND_ASSIGN(MatrixBase);
};

/// A class for storing matrices.
template<typename Real>
class Matrix : public MatrixBase<Real> {
 public:

  /// Empty constructor.
  Matrix();

  /// Basic constructor.
  Matrix(const MatrixIndexT r, const MatrixIndexT c,
         MatrixResizeType resize_type = kSetZero,
         MatrixStrideType stride_type = kDefaultStride):
      MatrixBase<Real>() { Resize(r, c, resize_type, stride_type); }

  /// Swaps the contents of *this and *other.  Shallow swap.
  void Swap(Matrix<Real> *other);

  /// Constructor from any MatrixBase. Can also copy with transpose.
  /// Allocates new memory.
  explicit Matrix(const MatrixBase<Real> & M,
                  MatrixTransposeType trans = kNoTrans);

  /// Same as above, but need to avoid default copy constructor.
  Matrix(const Matrix<Real> & M);  //  (cannot make explicit)

  /// Copy constructor: as above, but from another type.
  template<typename OtherReal>
  explicit Matrix(const MatrixBase<OtherReal> & M,
                    MatrixTransposeType trans = kNoTrans);

  /// read from stream.
  // Unlike one in base, allows resizing.
  void Read(std::istream & in, bool binary, bool add = false);

  /// Transpose the matrix.  Works for non-square
  /// matrices as well as square ones.
  void Transpose();

  /// Distructor to free matrices.
  ~Matrix() { Destroy(); }

  /// Sets matrix to a specified size (zero is OK as long as both r and c are
  /// zero).  The value of the new data depends on resize_type:
  ///   -if kSetZero, the new data will be zero
  ///   -if kUndefined, the new data will be undefined
  ///   -if kCopyData, the new data will be the same as the old data in any
  ///      shared positions, and zero elsewhere.
  ///
  /// You can set stride_type to kStrideEqualNumCols to force the stride
  /// to equal the number of columns; by default it is set so that the stride
  /// in bytes is a multiple of 16.
  ///
  /// This function takes time proportional to the number of data elements.
  void Resize(const MatrixIndexT r,
              const MatrixIndexT c,
              MatrixResizeType resize_type = kSetZero,
              MatrixStrideType stride_type = kDefaultStride);

  /// Assignment operator that takes MatrixBase.
  Matrix<Real> &operator = (const MatrixBase<Real> &other) {
    if (MatrixBase<Real>::NumRows() != other.NumRows() ||
        MatrixBase<Real>::NumCols() != other.NumCols())
      Resize(other.NumRows(), other.NumCols(), kUndefined);
    MatrixBase<Real>::CopyFromMat(other);
    return *this;
  }

  /// Assignment operator. Needed for inclusion in std::vector.
  Matrix<Real> &operator = (const Matrix<Real> &other) {
    if (MatrixBase<Real>::NumRows() != other.NumRows() ||
        MatrixBase<Real>::NumCols() != other.NumCols())
      Resize(other.NumRows(), other.NumCols(), kUndefined);
    MatrixBase<Real>::CopyFromMat(other);
    return *this;
  }


 private:
  /// Deallocates memory and sets to empty matrix (dimension 0, 0).
  void Destroy();

  /// Init assumes the current class contents are invalid (i.e. junk or have
  /// already been freed), and it sets the matrix to newly allocated memory with
  /// the specified number of rows and columns.  r == c == 0 is acceptable.  The data
  /// memory contents will be undefined.
  void Init(const MatrixIndexT r,
            const MatrixIndexT c,
            const MatrixStrideType stride_type);

};
/// @} end "addtogroup matrix_group"

/**
  Sub-matrix representation.
  Can work with sub-parts of a matrix using this class.
  Note that SubMatrix is not very const-correct-- it allows you to
  change the contents of a const Matrix.  Be careful!
*/

template<typename Real>
class SubMatrix : public MatrixBase<Real> {
 public:
  // Initialize a SubMatrix from part of a matrix; this is
  // a bit like A(b:c, d:e) in Matlab.
  // This initializer is against the proper semantics of "const", since
  // SubMatrix can change its contents.  It would be hard to implement
  // a "const-safe" version of this class.
  SubMatrix(const MatrixBase<Real>& T,
            const MatrixIndexT ro,  // row offset, 0 < ro < NumRows()
            const MatrixIndexT r,   // number of rows, r > 0
            const MatrixIndexT co,  // column offset, 0 < co < NumCols()
            const MatrixIndexT c);   // number of columns, c > 0

  // This initializer is mostly intended for use in CuMatrix and related
  // classes.  Be careful!
  SubMatrix(Real *data,
            MatrixIndexT num_rows,
            MatrixIndexT num_cols,
            MatrixIndexT stride);

  ~SubMatrix<Real>() {}

  /// This type of constructor is needed for Range() to work [in Matrix base
  /// class]. Cannot make it explicit.
  SubMatrix<Real> (const SubMatrix &other):
  MatrixBase<Real> (other.data_, other.num_cols_, other.num_rows_,
                    other.stride_) {}

 private:
  /// Disallow assignment.
  SubMatrix<Real> &operator = (const SubMatrix<Real> &other);
};
/// @} End of "addtogroup matrix_funcs_io".

/// \addtogroup matrix_funcs_scalar
/// @{

// Some declarations.  These are traces of products.


template<typename Real>
inline void AssertEqual(const MatrixBase<Real> &A, const MatrixBase<Real> &B,
                        float tol = 0.01) {
  KALDI_ASSERT(A.ApproxEqual(B, tol));
}

/// @} end "addtogroup matrix_funcs_scalar"


/// \addtogroup matrix_funcs_io
/// @{
template<typename Real>
std::ostream & operator << (std::ostream & Out, const MatrixBase<Real> & M);

template<typename Real>
std::istream & operator >> (std::istream & In, MatrixBase<Real> & M);

// The Matrix read allows resizing, so we override the MatrixBase one.
template<typename Real>
std::istream & operator >> (std::istream & In, Matrix<Real> & M);


template<typename Real>
bool SameDim(const MatrixBase<Real> &M, const MatrixBase<Real> &N) {
  return (M.NumRows() == N.NumRows() && M.NumCols() == N.NumCols());
}
/// @} end of \addtogroup matrix_funcs_io

}  // namespace kaldi



// we need to include the implementation and some
// template specializations.
#include "matrix/kaldi-matrix-inl.h"


#endif  // KALDI_MATRIX_KALDI_MATRIX_H_
