The instructions of src of the Kaldi-based standalone pitch extractor toolkit.

The core idea
=================================================

Our algorithm is a highly modified version of the getf0 algorithm. It is based
on finding lag values that maximize the Normalized Cross Correlation Function.
Probably the most important change from getf0 is that rather than making hard
decisions about voicing on each frame, we treat all frames as voiced and allow
the Viterbi search to naturally interpolate across unviced regions.


The structure of src directory
=================================================

"featbin/" There are five key binaries in this directory. They are the whole tools
in the toolkit. As a user, you need to select one tool which statifies your requirement
and then conduct it. We will introduce the five tools separately in below.

"feat/" It provide the necessary "class" and "function" for extracting pitch feature.

"matrix/" The basic matrix or vector operations are implemented in this directory.

Other directories are used to supply the toolkit.


The specific introduction of tools
=================================================

The five tools in featbin/ directory are "compute-and-process-kaldi-pitch-feats",
"compute-kaldi-pitch-feats", "process-kaldi-pitch-feats", "interpolate-pitch"
and "process-pitch-feats."
They can be divided into two parts. The first part contains the first three tools,
and the remainings belong to the second part.

1. The first part
   "compute-and-process-kaldi-pitch-feats" = "compute-kaldi-pitch-feats" +
                                             "process-kaldi-pitch-feats"
   
   "compute-and-process-kaldi-pitch-feats" applies pitch extractor and pitch
   post-processor, starting from wav input. The output information is as follows. 
   
   At most, the output contains a 4-tuple. It is (POV-feature, 
   Normalized-log-pitch, Delta-pitch, Raw-log-pitch). You can select from four 
   features by setting the boolean options (--add-pov-feature,
   --add-normalized-log-pitch, --add-delta-pitch and --add-raw-log-pitch).
   
   Default setup produces 3-tuple features which consists of (POV-feature,
   pitch-feature, delta-pitch-feautre).
   POV-feautre is warped NCCF.
   Pitch-feature is log-pitch with POV-weighted mean subtraction over 1.5s window.
   Delta-pitch-feature is delta feature computed on raw log pitch.


   "compute-kaldi-pitch-feats" applies pitch extractor, starting from wav input.
   Output is 2-tuple features consisting of (NCCF, pitch).
   Nccf is between -1 and 1, and higher for voiced frames.

   "process-kaldi-pitch-feats" applies post-process pitch features, consisting
   of pitch and NCCF, into features suitable for input to ASR system. See above
   to get the information about output.

2. The second part
   "Interpolate-pitch" is a rather special-purpose program which processes 
   2-dimensional features consisting of (prob-of-voicing, pitch).  By default 
   we do model-based pitch smoothing and interpolation.

   "process-pitch-feats" is a rather special-purpose program which processes 
   2-dimensional features consisting of (prob-of-voicing, pitch) into something
   suitable to put into a speech recognizer.


The main handset parameters
=================================================

Like most pitch extraction algorithms, our algorithm has a number of parameters
that were set by hand. The following table shows the information

-------------------------------------------------------------------------------
Parameter             Value       Explanation
-------------------------------------------------------------------------------
min-f0                50          Minimum possible frequency value (Hz)
max-f0                400         Maximum possible frequency value (Hz)
window-width          0.025       Length in seconds of window used for NCCF
window-shift          0.01        Frame-shift, in seconds (should match
                                  that used for baseline features e.g. PLP)
soft-min-f0           10          Minimum f0, applied in soft way;
                                  must not exceed min-f0.
nccf-ballast          0.625       Increasing this factor reduces NCCF for quiet
                                  frames, helping ensure pitch continuity
                                  in unvoiced regions
penalty-factor        0.1         Factor that penalizes frequency change
delta-pitch           0.005       Smallest relative change in pitch
                                  that our algorithm measures
lowpass-cutoff        1000        Low-pass cutoff that we apply to
                                  the raw signal
lowpass-filter-width  2           Integer that determines filter width
                                  of low-pass filter (more gives wider filter with
                                  sharper cutoff)
resample-frequency    4000        Sample frequency for NCCF;
                                  must exceed twice lowpass-cutoff.
upsample-filter-width 5           Integer that determines filter width
                                  when upsampling NCCF
-------------------------------------------------------------------------------
