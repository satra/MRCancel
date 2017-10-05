function yout = soundframer(y,fs,framelen,evalfunc,dispfunc,params)
% SOUNDFRAMER Evaluates a process on framed signals

persistent hproc

if nargin<6,
    params = struct([]);
end

Ny = numel(y);
yout = zeros(size(y));
plot_bar = 1;
for f0=1:framelen:length(y),
    frameidx = f0:min(f0+framelen-1,length(y));
    frame = y(frameidx);
    if (f0 == 1)
        if (plot_bar)
        hproc = uiwaitbar('Processing');
        end
        feval(evalfunc,'init',frame,fs,f0,params);
    end
    [outframe] = feval(evalfunc,'process',frame,fs,f0,params);
    if mod(f0-1,80*framelen)==0 && plot_bar,
        uiwaitbar(frameidx(end)/Ny,hproc);
    end
    yout(frameidx,1) = outframe;
%     plot([y(1:frameidx(end)),yout(params.delay+[1:frameidx(end)])]);
%     drawnow;
    if ~isempty(dispfunc),
        feval(dispfunc,yout,frameidx,f0);
    end
end
if (plot_bar)
   delete(hproc);
end