function clean_recordings2(in_file, out_file, TR, slices, rmsupdate, channel, onset, windowext, rmsthresh)
% function clean_recodings(in_file, out_file, TR, slices, rmsupdate, channel, onset, windowext, rmsthresh)
%
% Removes scanner noise from voice recordings
%    clean_recordings('path_to_input', 'path_to_output', TR, slices, rmsupdate=0.007, channel=1, onset=0, windowext=20, rmsthresh=0.5);
%
% satra@mit.edu


% Associated recall experiment

if nargin<5, rmsupdate = 0.007; end
if nargin<6, channel = 1; end
if nargin<7, onset = 0; end
if nargin<8, windowext = 20; end
if nargin<5, rmsthresh = 0.5; end
 
global templates rmsres_store cmax_store template
templates = [];
template = [];
rmsres_store = [];
cmax_store = [];


% Read the file
[y,fs] = audioread(in_file);
y = y((onset*fs):end, channel); %(100*fs), 2); % zeros(4*fs, 1)];
framelen = round(0.025*fs); %round(0.025*fs);

%% Parameters
p.TR = TR; % Time of repetition
p.slices = slices;

p.window = p.TR/p.slices;             %approximate template window length [ms]
p.windowext = round(windowext);            %distortion in window length [ms]
p.template_corr1 = 0.997; %0.997;     %correlation threshold for noise template match
p.template_corr2 = 0.97;     %correlation threshold when noise template
%selected as a result of old correlation>new correlation
p.rmsthresh  = rmsthresh;         %rms threshold on the normalized signal for
%either looking for signal or taking any
%corrective action
p.corr_remove = 0.9; %0.55        %remove noise only when old correlation is greater than this value
p.nc_param = 1;           %weight on spectral noise subtraction [1 == exact subtraction]
%p.wt_type  = 1;             %shape of weighting
%0,1 : set wt to these values
%2: phase signal_transform*conj(noise_transform)
%3: funky wt using unwrap
%4: linear asymmetric wt
p.cutoff = 4000;        %digital filter cutoff [Hz]
%p.rmsupdate = 0.005; %0.025;         %template updated when rms value of residual is less than this value
p.alpha = 0.5;              %template = p.alpha*old_template +(1-p.alpha)*new_template;
p.delay = floor(framelen+(2*(p.window+p.windowext))*1e-3*fs); %computed delay

%         p.wt_type = 1;p.rmsupdate = 0.0003; %0.025;         %template updated when rms value of residual is less than this value
%05   p.wt_type = 6;p.rmsupdate = 0.008; %0.025;         %template updated when rms value of residual is less than this value
%04     p.wt_type = 6;p.rmsupdate = 0.021; %0.025;         %template updated when rms value of residual is less than this value
%06
p.wt_type = 6;p.rmsupdate = rmsupdate; %0.025;         %template updated when rms value of residual is less than this value

%% setup signal
maxy = max(abs(y));
y = y/maxy;
y = y(:);

%% reduce noise and write output
yout = soundframer(y,fs,framelen,@rn10,'',p);
audiowrite(out_file, yout.*(yout<4), fs);

