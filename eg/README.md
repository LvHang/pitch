The instructions of eg of the Kaldi-based standalone pitch extractor toolkit.

A simple example
=================================================

In the eg/ directory, we provide a simple example for using the toolkit. See
run.sh for details.

The tools are based on Kaldi I/O mechanisms. For the rspecifier or wspecifier,
you need to use "ark" or "scp" symbol to assign the input files and output files.
Please check Kaldi website for more information about Kaldi I/O mechanisms.
(http://kaldi-asr.org/doc/io.html)

When you use "compute-and-process-kaldi-pitch-feats" or "compute-kaldi-pitch-feats",
you need to prepare the "wav.scp" file as input "<feats-rspecifier>", whose format
is <recording-id> <extended-filename>, where the "extended-filename" may be an 
actual filename, or as in this case, a command that extracts a wav-format file. 
The pipe symbol on the end of the extended-filename specifies that it is to be 
interpreted as a pipe.

When you use "process-kaldi-pitch-feats", "interpolate-pitch" or "process-pitch-feats".
The input is "<feats-rspecifier>", which is the output file of 
"compute-and-process-kaldi-pitch-feats" or "compute-kaldi-pitch-feats". The 
"<feats-rspecifier>" should be "scp" file or "ark" file.
The format of "scp" file is <recording-id> <information-position-in-ark-file>.
It is treated as a index form. The format of "ark" file is <recording-id> <infromation>,
which stores the specific data.


The specific introduction of tools
=================================================

There are five tools in this toolkit, which are named 
"compute-and-process-kaldi-pitch-feats", "compute-kaldi-pitch-feats", 
"process-kaldi-pitch-feats", "interpolate-pitch" and "process-pitch-feats."

Please check `src/README.md` for the informaiton of the five tools. For the 
information of options of each tool, please check the usage of a specific tool.
