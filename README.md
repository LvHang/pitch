The Kaldi-based Standalone Pitch Extractor Toolkit

OverView
=================================================

1. The standalone pitch extractor toolkit is based on Kaldi toolkit
   (https://github.com/kaldi-asr/kaldi)

2. The implementation of the standalone pitch extractor toolkit is based on
   the paper--A Pitch Extraction Algorithm Tuned For Automatic Speech Recognition,
   whose author is Pegah Ghahremani, Daniel Povey and etc, in technology.
   (http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6854049)

   Our algorithm is a highly modified version of the getf0 algorithm. It is based
   on finding lag values that maximize the Normalized Cross Correlation Function.
   Probably the most important change from getf0 is that rather than making hard
   decisions about voicing on each frame, we treat all frames as voiced and allow
   the Viterbi search to naturally interpolate across unviced regions.

3. The toolkit does not depend on any other toolkits, so all the source codes
   are located in src/. For more information about source codes, see `src/README.md`

   To build the toolkit: see `src/INSTALL`. These instructions are valid for UNIX
   systems including various flavors of Linux; Darwin; and Cygwin (has not been 
   tested on more "exotic" varieties of UNIX).

   To run the example system bilds, see `eg/README.md`

   If you encounter problems (and you probably will), please do not hesitate to
   contact the developers (see below). In addition to specific qestions, please 
   let us know if there are specific aspects of the project that you feel could be
   improved, that you find confusing, etc., and which missing features you most
   wish it had.


Development pattern for contributors
=================================================

1. [Create a personal fork](https://help.github.com/articles/fork-a-repo/)
   of the [main Kaldi repository](https://github.com/kaldi-asr/kaldi) in GitHub.
2. Make your changes in a named branch different from `master`, e.g. you create
   a branch `my-awesome-feature`.
3. [Generate a pull request](https://help.github.com/articles/creating-a-pull-request/)
   through the Web interface of GitHub.
4. As a general rule, please follow [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
   There are a [few exceptions in Kaldi](http://kaldi-asr.org/doc/style.html).
   You can use the [Google's cpplint.py](https://raw.githubusercontent.com/google/styleguide/gh-pages/cpplint/cpplint.py)
   to verify that your code is free of basic mistakes.


Authors
=================================================
Daniel Povey (dpovey@gmail.com)
Pegah Ghahremani (pghahre1@jhu.edu)
Hang Lyu (hanglv@nwpu-aslp.org)


Platform specific notes
=================================================

### PowerPC 64bits little-endian (ppc64le)

- The tool is expected to work out of the box in RHEL >= 7 and Ubuntu >= 16.04

### Android
### (TODO: Test)

- The tool supports cross compiling for Android using Android NDK, clang++
- See [this blog post](http://jcsilva.github.io/2017/03/18/compile-kaldi-android/)
  for details.
