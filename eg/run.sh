#!/bin/bash

# Begin configuration section.
nj=4
cmd=utils/run.pl
pitch_config=conf/pitch.conf
# End configuration section.

echo "$0 $@"  # Print the command line for logging

if [ -f path.sh ]; then . ./path.sh; fi
. utils/parse_options.sh || exit 1;

if [ $# -lt 2 ] || [ $# -gt 4 ]; then
   echo "Usage: $0 [options] <tool-index> <data-dir> [<log-dir> [<output-dir>] ]";
   echo "e.g.: $0 data/ data/log data/information"
   echo "Note: <log-dir> defaults to <data-dir>/log, and <output-dir> defaults to <data-dir>/information"
   echo "Please modify the corresponding command or use config file, when you want to use non-default options."
   echo "Options: "
   echo "  tool-index <1|2|3|4|5>             # 1--compute-and-proess-kaldi-pitch-feats"
   echo "                                     # 2--compute-kaldi-pitch-feats"
   echo "                                     # 3--process-kaldi-pitch-feats"
   echo "                                     # 4--interpolate-pitch"
   echo "                                     # 5--process-pitch-feats"
   echo "  --pitch-config <pitch-config-file> # config passed to compute-and-process-kaldi-pitch-feats etc. "
   echo "                                     # default is conf/pitch.conf"
   echo "   "
   exit 1;
fi

index=$1
data=$2
if [ $# -ge 3 ]; then
  logdir=$3
else
  logdir=$data/log
fi
if [ $# -ge 4 ]; then
  pitch_dir=$4
else
  pitch_dir=$data/information
fi

# make $pitch_dir an absolute pathname.
pitch_dir=`perl -e '($dir,$pwd)= @ARGV; if($dir!~m:^/:) { $dir = "$pwd/$dir"; } print $dir; ' $pitch_dir ${PWD}`

# use "name" as part of name of the archive.
name=`basename $data`

mkdir -p $pitch_dir || exit 1;
mkdir -p $logdir || exit 1;

if [[ $index -lt 1 || $index -gt 5 ]]; then
   echo "Please select a tool index from 1 to 5."
   exit 1
fi

if [ $index == 1 ]; then
  scp=$data/wav.scp
  if [ ! -f $scp ]; then
    echo "compute-and-process-kaldi-pitch-feats needs the input wav.scp"
    exit 1;
  fi
  $cmd $logdir/compute-and-process-kaldi-pitch-feats.log \
    compute-and-process-kaldi-pitch-feats scp:$scp \
    ark,scp:$pitch_dir/feats.ark,$pitch_dir/feats.scp || exit 1
  cp $pitch_dir/feats.scp $data/feats.scp
fi
if [ $index == 2 ]; then
  scp=$data/wav.scp
  if [ ! -f $scp ]; then
    echo "compute-kaldi-pitch-feats needs the input wav.scp"
    exit 1;
  fi
  $cmd $logdir/compute-kaldi-pitch-feats.log \
    compute-kaldi-pitch-feats scp:$scp \
    ark,scp:$pitch_dir/feats.ark,$pitch_dir/feats.scp || exit 1
  cp $pitch_dir/feats.scp $data/feats.scp
fi
if [ $index == 3 ]; then
  scp=$data/wav.scp
  if [ ! -f $data/wav.scp ]; then
    echo "process-kaldi-pitch-feats needs the input wav.scp"
    exit 1;
  fi
  $cmd $logdir/process-kaldi-pitch-feats.log \
    process-kaldi-pitch-feats ark:"compute-kaldi-pitch-feats scp:$scp ark:- |" \
    ark,scp:$pitch_dir/feats.ark,$pitch_dir/feats.scp || exit 1
  cp $pitch_dir/feats.scp $data/feats.scp
fi
if [ $index == 4 ]; then
  scp=$data/wav.scp
  if [ ! -f $data/wav.scp ]; then
    echo "interpolate-pitch needs the input wav.scp"
    exit 1;
  fi
  $cmd $logdir/interpolate-pitch.log \
    interpolate-pitch ark:"compute-kaldi-pitch-feats scp:$scp ark:- |" \
    ark,scp:$pitch_dir/feats.ark,$pitch_dir/feats.scp || exit 1
  cp $pitch_dir/feats.scp $data/feats.scp
fi
if [ $index == 5 ]; then
  scp=$data/wav.scp
  if [ ! -f $data/wav.scp ]; then
    echo "process-pitch-feats needs the input wav.scp"
    exit 1;
  fi
  $cmd $logdir/process-pitch-feats.log \
    process-pitch-feats \
    ark:"compute-kaldi-pitch-feats scp:$scp ark:- |" \
    ark,scp:$pitch_dir/feats.ark,$pitch_dir/feats.scp || exit 1
  cp $pitch_dir/feats.scp $data/feats.scp
fi 

echo "Succeeded"
