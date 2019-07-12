The intent of this script is to remove MR scanner pulse from audio recordings
of a participant speaking in the scanner during EPI acquisition.

For details and citing this software see: https://arxiv.org/abs/1207.5827v1

This program has been tested to work with MATLAB R2018a.

### Compile the mex files.
- mex testmult.c
- mex create_template03.c


### Usage:

clean_recordings('path_to_input', 'path_to_output', TR, slices, rmsupdate=0.007, channel=1, onset=0, windowext=20, rmsthresh=0.5);

TR: Repetition Time in milliseconds
slices: Total number of slices or slice groups acquired  
rmsupdate: Threshold below which a template is updated. Start with the default and slowly increase it if templates are not being updated.
channel: which channel to use in a stereo recording
onset: onset of EPI pulses
windowext: slop factor on template length (TR/slices +/- windowext)
rmsthresh: Threshold to use to ensure there is sound in the buffer

Example: 
clean_recordings2('in_filename.wav', 'output.wav', 2000, 32, 0.03, 2, 35.5, 5, 0.2);

If you have questions please post in the issues of this repository.