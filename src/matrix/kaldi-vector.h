// matrix/kaldi-vector.h

// Copyright 2009-2012   Ondrej Glembek;  Microsoft Corporation;  Lukas Burget;
//                       Saarland University (Author: Arnab Ghoshal);
//                       Ariya Rastrow;  Petr Schwarz;  Yanmin Qian;
//                       Karel Vesely;  Go Vivace Inc.;  Arnab Ghoshal
//                       Wei Shi;
//                2015   Guoguo Chen

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

#ifndef KALDI_MATRIX_KALDI_VECTOR_H_
#define KALDI_MATRIX_KALDI_VECTOR_H_ 1

#include "matrix/matrix-common.h"

namespace kaldi {

/// \addtogroup matrix_group
/// @{

///  Provides a vector abstraction class.
///  This class provides a way to work with vectors in kaldi.
///  It encapsulates basic operations and memory optimizations.
template<typename Real>
class VectorBase {
 public:
  /// Set vector to all zeros.
  void SetZero();

  /// Returns true if matrix is all zeros.
  bool IsZero(Real cutoff = 1.0e-06) const;     // replace magic number

  /// Set all members of a vector to a specified value.
  void Set(Real f);

  /// Set vector to random normally-distributed noise.
  void SetRandn();

  /// Returns the  dimension of the vector.
  inline MatrixIndexT Dim() const { return dim_; }

  /// Returns the size in memory of the vector, in bytes.
  inline MatrixIndexT SizeInBytes() const { return (dim_*sizeof(Real)); }

  /// Returns a pointer to the start of the vector's data.
  inline Real* Data() { return data_; }

  /// Returns a pointer to the start of the vector's data (const).
  inline const Real* Data() const { return data_; }

  /// Indexing  operator (const).
  inline Real operator() (MatrixIndexT i) const {
    KALDI_PARANOID_ASSERT(static_cast<UnsignedMatrixIndexT>(i) <
                 static_cast<UnsignedMatrixIndexT>(dim_));
    return *(data_ + i);
  }

  /// Indexing operator (non-const).
  inline Real & operator() (MatrixIndexT i) {
    KALDI_PARANOID_ASSERT(static_cast<UnsignedMatrixIndexT>(i) <
                 static_cast<UnsignedMatrixIndexT>(dim_));
    return *(data_ + i);
  }

  /** @brief Returns a sub-vector of a vector (a range of elements).
   *  @param o [in] Origin, 0 < o < Dim()
   *  @param l [in] Length 0 < l < Dim()-o
   *  @return A SubVector object that aliases the data of the Vector object.
   *  See @c SubVector class for details   */
  SubVector<Real> Range(const MatrixIndexT o, const MatrixIndexT l) {
    return SubVector<Real>(*this, o, l);
  }

  /** @brief Returns a const sub-vector of a vector (a range of elements).
   *  @param o [in] Origin, 0 < o < Dim()
   *  @param l [in] Length 0 < l < Dim()-o
   *  @return A SubVector object that aliases the data of the Vector object.
   *  See @c SubVector class for details   */
  const SubVector<Real> Range(const MatrixIndexT o,
                              const MatrixIndexT l) const {
    return SubVector<Real>(*this, o, l);
  }

  /// Copy data from another vector (must match own size).
  void CopyFromVec(const VectorBase<Real> &v);

  /// Copy data from another vector of different type (double vs. float)
  template<typename OtherReal>
  void CopyFromVec(const VectorBase<OtherReal> &v);

  /// Applies floor to all elements. Returns number of elements floored.
  MatrixIndexT ApplyFloor(Real floor_val);

  /// Applies floor to all elements. Returns number of elements floored.
  MatrixIndexT ApplyFloor(const VectorBase<Real> &floor_vec);

  /// Take all  elements of vector to a power.
  void ApplyPow(Real power);

  /// Compute the p-th norm of the vector.
  Real Norm(Real p) const;

  /// Returns true if ((*this)-other).Norm(2.0) <= tol * (*this).Norm(2.0).
  bool ApproxEqual(const VectorBase<Real> &other, float tol = 0.01) const;

  /// Invert all elements.
  void InvertElements();

  /// Add vector : *this = *this + alpha * rv (with casting between floats and
  /// doubles)
  template<typename OtherReal>
  void AddVec(const Real alpha, const VectorBase<OtherReal> &v);

  /// Add vector : *this = *this + alpha * rv^2  [element-wise squaring].
  void AddVec2(const Real alpha, const VectorBase<Real> &v);

  /// Add vector : *this = *this + alpha * rv^2  [element-wise squaring],
  /// with casting between floats and doubles.
  template<typename OtherReal>
  void AddVec2(const Real alpha, const VectorBase<OtherReal> &v);

  /// Add matrix times vector : this <-- beta*this + alpha*M*v.
  void AddMatVec(const Real alpha, const MatrixBase<Real> &M,
                 const MatrixTransposeType trans,  const VectorBase<Real> &v,
                 const Real beta); // **beta previously defaulted to 0.0**

  /// Multipy element-by-element by another vector.
  void MulElements(const VectorBase<Real> &v);
  /// Multipy element-by-element by another vector of different type.
  template<typename OtherReal>
  void MulElements(const VectorBase<OtherReal> &v);

  /// Add a constant to each element of a vector.
  void Add(Real c);

  /// Add element-by-element product of vectlrs:
  //  this <-- alpha * v .* r + beta*this .
  void AddVecVec(Real alpha, const VectorBase<Real> &v,
                 const VectorBase<Real> &r, Real beta);

  /// Multiplies all elements by this constant.
  void Scale(Real alpha);

  /// Returns the maximum value of any element, or -infinity for the empty vector.
  Real Max() const;

  /// Returns the maximum value of any element, and the associated index.
  /// Error if vector is empty.
  Real Max(MatrixIndexT *index) const;

  /// Returns the minimum value of any element, or +infinity for the empty vector.
  Real Min() const;

  /// Returns the minimum value of any element, and the associated index.
  /// Error if vector is empty.
  Real Min(MatrixIndexT *index) const;

  /// Returns sum of the elements
  Real Sum() const;

  /// Does *this = alpha * (sum of rows of M) + beta * *this.
  void AddRowSumMat(Real alpha, const MatrixBase<Real> &M, Real beta = 1.0);

  /// Add the diagonal of a matrix times itself:
  /// *this = diag(M M^T) +  beta * *this (if trans == kNoTrans), or
  /// *this = diag(M^T M) +  beta * *this (if trans == kTrans).
  void AddDiagMat2(Real alpha, const MatrixBase<Real> &M,
                   MatrixTransposeType trans = kNoTrans, Real beta = 1.0);

  /// Reads from C++ stream (option to add to existing contents).
  /// Throws exception on failure
  void Read(std::istream & in, bool binary, bool add = false);

  /// Writes to C++ stream (option to write in binary).
  void Write(std::ostream &Out, bool binary) const;

  friend class VectorBase<double>;
  friend class VectorBase<float>;
 protected:
  /// Destructor;  does not deallocate memory, this is handled by child classes.
  /// This destructor is protected so this object can only be deleted via a child.
  ~VectorBase() {}

  /// Empty initializer, corresponds to vector of zero size.
  explicit VectorBase(): data_(NULL), dim_(0) {
    KALDI_ASSERT_IS_FLOATING_TYPE(Real);
  }

  /// data memory area
  Real* data_;
  /// dimension of vector
  MatrixIndexT dim_;
  KALDI_DISALLOW_COPY_AND_ASSIGN(VectorBase);
}; // class VectorBase

/** @brief A class representing a vector.
 *
 *  This class provides a way to work with vectors in kaldi.
 *  It encapsulates basic operations and memory optimizations.  */
template<typename Real>
class Vector: public VectorBase<Real> {
 public:
  /// Constructor that takes no arguments.  Initializes to empty.
  Vector(): VectorBase<Real>() {}

  /// Constructor with specific size.  Sets to all-zero by default
  /// if set_zero == false, memory contents are undefined.
  explicit Vector(const MatrixIndexT s,
                  MatrixResizeType resize_type = kSetZero)
      : VectorBase<Real>() {  Resize(s, resize_type);  }

  /// Copy constructor.  The need for this is controversial.
  Vector(const Vector<Real> &v) : VectorBase<Real>()  { //  (cannot be explicit)
    Resize(v.Dim(), kUndefined);
    this->CopyFromVec(v);
  }

  /// Copy-constructor from base-class, needed to copy from SubVector.
  explicit Vector(const VectorBase<Real> &v) : VectorBase<Real>() {
    Resize(v.Dim(), kUndefined);
    this->CopyFromVec(v);
  }

  /// Type conversion constructor.
  template<typename OtherReal>
  explicit Vector(const VectorBase<OtherReal> &v): VectorBase<Real>() {
    Resize(v.Dim(), kUndefined);
    this->CopyFromVec(v);
  }

  /// Swaps the contents of *this and *other.  Shallow swap.
  void Swap(Vector<Real> *other);

  /// Destructor.  Deallocates memory.
  ~Vector() { Destroy(); }

  /// Read function using C++ streams.  Can also add to existing contents
  /// of matrix.
  void Read(std::istream & in, bool binary, bool add = false);

  /// Set vector to a specified size (can be zero).
  /// The value of the new data depends on resize_type:
  ///   -if kSetZero, the new data will be zero
  ///   -if kUndefined, the new data will be undefined
  ///   -if kCopyData, the new data will be the same as the old data in any
  ///      shared positions, and zero elsewhere.
  /// This function takes time proportional to the number of data elements.
  void Resize(const MatrixIndexT dim, MatrixResizeType resize_type = kSetZero);

  /// Assignment operator, protected so it can only be used by std::vector
  Vector<Real> &operator = (const Vector<Real> &other) {
    Resize(other.Dim(), kUndefined);
    this->CopyFromVec(other);
    return *this;
  }

  /// Assignment operator that takes VectorBase.
  Vector<Real> &operator = (const VectorBase<Real> &other) {
    Resize(other.Dim(), kUndefined);
    this->CopyFromVec(other);
    return *this;
  }
 private:
  /// Init assumes the current contents of the class are invalid (i.e. junk or
  /// has already been freed), and it sets the vector to newly allocated memory
  /// with the specified dimension.  dim == 0 is acceptable.  The memory contents
  /// pointed to by data_ will be undefined.
  void Init(const MatrixIndexT dim);

  /// Destroy function, called internally.
  void Destroy();

};


/// Represents a non-allocating general vector which can be defined
/// as a sub-vector of higher-level vector [or as the row of a matrix].
template<typename Real>
class SubVector : public VectorBase<Real> {
 public:
  /// Constructor from a Vector or SubVector.
  /// SubVectors are not const-safe and it's very hard to make them
  /// so for now we just give up.  This function contains const_cast.
  SubVector(const VectorBase<Real> &t, const MatrixIndexT origin,
            const MatrixIndexT length) : VectorBase<Real>() {
    // following assert equiv to origin>=0 && length>=0 &&
    // origin+length <= rt.dim_
    KALDI_ASSERT(static_cast<UnsignedMatrixIndexT>(origin)+
                 static_cast<UnsignedMatrixIndexT>(length) <=
                 static_cast<UnsignedMatrixIndexT>(t.Dim()));
    VectorBase<Real>::data_ = const_cast<Real*> (t.Data()+origin);
    VectorBase<Real>::dim_   = length;
  }

  /// Copy constructor
  SubVector(const SubVector &other) : VectorBase<Real> () {
    // this copy constructor needed for Range() to work in base class.
    VectorBase<Real>::data_ = other.data_;
    VectorBase<Real>::dim_ = other.dim_;
  }

  /// Constructor from a pointer to memory and a length.  Keeps a pointer
  /// to the data but does not take ownership (will never delete).
  SubVector(Real *data, MatrixIndexT length) : VectorBase<Real> () {
    VectorBase<Real>::data_ = data;
    VectorBase<Real>::dim_   = length;
  }


  /// This operation does not preserve const-ness, so be careful.
  SubVector(const MatrixBase<Real> &matrix, MatrixIndexT row) {
    VectorBase<Real>::data_ = const_cast<Real*>(matrix.RowData(row));
    VectorBase<Real>::dim_   = matrix.NumCols();
  }

  ~SubVector() {}  ///< Destructor (does nothing; no pointers are owned here).

 private:
  /// Disallow assignment operator.
  SubVector & operator = (const SubVector &other) {}
};

/// @} end of "addtogroup matrix_group"
/// \addtogroup matrix_funcs_io
/// @{
/// Output to a C++ stream.  Non-binary by default (use Write for
/// binary output).
template<typename Real>
std::ostream & operator << (std::ostream & out, const VectorBase<Real> & v);

/// Input from a C++ stream.  Will automatically read text or
/// binary data from the stream.
template<typename Real>
std::istream & operator >> (std::istream & in, VectorBase<Real> & v);

/// Input from a C++ stream. Will automatically read text or
/// binary data from the stream.
template<typename Real>
std::istream & operator >> (std::istream & in, Vector<Real> & v);
/// @} end of \addtogroup matrix_funcs_io

/// \addtogroup matrix_funcs_scalar
/// @{


template<typename Real>
bool ApproxEqual(const VectorBase<Real> &a,
                 const VectorBase<Real> &b, Real tol = 0.01) {
  return a.ApproxEqual(b, tol);
}

template<typename Real>
inline void AssertEqual(VectorBase<Real> &a, VectorBase<Real> &b,
                        float tol = 0.01) {
  KALDI_ASSERT(a.ApproxEqual(b, tol));
}


/// Returns dot product between v1 and v2.
template<typename Real>
Real VecVec(const VectorBase<Real> &v1, const VectorBase<Real> &v2);

template<typename Real, typename OtherReal>
Real VecVec(const VectorBase<Real> &v1, const VectorBase<OtherReal> &v2);


/// @} End of "addtogroup matrix_funcs_scalar"


}  // namespace kaldi

// we need to include the implementation
#include "matrix/kaldi-vector-inl.h"



#endif  // KALDI_MATRIX_KALDI_VECTOR_H_

