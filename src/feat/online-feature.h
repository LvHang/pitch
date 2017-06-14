// feat/online-feature.h

// Copyright 2013   Johns Hopkins University (author: Daniel Povey)
//           2014   Yanqing Sun, Junjie Wang,
//                  Daniel Povey, Korbinian Riedhammer

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


#ifndef KALDI_FEAT_ONLINE_FEATURE_H_
#define KALDI_FEAT_ONLINE_FEATURE_H_

#include <string>
#include <vector>
#include <deque>

#include "matrix/matrix-lib.h"
#include "util/common-utils.h"
#include "base/kaldi-error.h"
#include "feat/feature-functions.h"
#include "itf/online-feature-itf.h"

namespace kaldi {
/// @addtogroup  onlinefeat OnlineFeatureExtraction
/// @{


/// This class takes a Matrix<BaseFloat> and wraps it as an
/// OnlineFeatureInterface: this can be useful where some earlier stage of
/// feature processing has been done offline but you want to use part of the
/// online pipeline.
class OnlineMatrixFeature: public OnlineFeatureInterface {
 public:
  /// Caution: this class maintains the const reference from the constructor, so
  /// don't let it go out of scope while this object exists.
  explicit OnlineMatrixFeature(const MatrixBase<BaseFloat> &mat): mat_(mat) { }

  virtual int32 Dim() const { return mat_.NumCols(); }

  virtual BaseFloat FrameShiftInSeconds() const {
    return 0.01f;
  }

  virtual int32 NumFramesReady() const { return mat_.NumRows(); }

  virtual void GetFrame(int32 frame, VectorBase<BaseFloat> *feat) {
    feat->CopyFromVec(mat_.Row(frame));
  }

  virtual bool IsLastFrame(int32 frame) const {
    return (frame + 1 == mat_.NumRows());
  }


 private:
  const MatrixBase<BaseFloat> &mat_;
};

/// @} End of "addtogroup onlinefeat"
}  // namespace kaldi

#endif  // KALDI_FEAT_ONLINE_FEATURE_H_
