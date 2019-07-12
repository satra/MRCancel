function [outframe] = rn10(action,frame,fs,framecount,p)
% NOTES:
% RMS needs to be callibrated
%

% create persistent buffer
global template rmsres_store cmax_store templates
persistent delay buffer buffer2 endmark maxc midx maxc2 tidx spidx
%nextindex

plot_debug = 0;

N = length(frame);
switch action,
    case 'init'
        if nargin<5  | isempty(p)
            %% Parameters
            p.window = 100;             %approximate template window length [ms]
            p.windowext = 4;            %distortion in window length [ms]
            p.template_corr1 = 0.99;     %correlation threshold for noise template match
            p.template_corr2 = 0.95;     %correlation threshold when noise template
            %selected as a result of old correlation>new correlation
            p.rmsthresh  = 0.1;         %rms threshold on the normalized signal for
            %either looking for signal or taking any
            %corrective action
            p.corr_remove = 0.55;        %remove noise only when old correlation is greater than this value
            p.nc_param = 1.2;           %weight on spectral noise subtraction [1 == exact subtraction]
            p.wt_type  = 1;             %shape of weighting
            %0,1 : set wt to these values
            %2: phase signal_transform*conj(noise_transform)
            %3: funky wt using unwrap
            %4: linear asymmetric wt
            p.cutoff = 0; %5000;        %digital filter cutoff [Hz]
            p.rmsupdate = 0.02;         %template updated when rms value of residual is less than this value
            p.alpha = 0.2;              %weighting of current template added
            p.delay = floor(framelen+(2*(p.window+p.windowext))*1e-3*fs); %computed delay
            fprintf('using default parameters\n');
        else
            fprintf('Not using default parameters\n');
        end
        delay = p.delay;
        buffer = zeros(delay,1);
        buffer2 = zeros(delay,1);
        endmark = inf;
        maxc = 2;
        maxc2 = 2;
        midx = 1;
        template = [];
        outframe = [];
        clear c;
        clear testmult;
        clear create_template03;
        return;
    case 'process',
        outframe = buffer(1:N);

        % Add to buffer
        buffer = push_q(buffer,frame,N);
        buffer2 = push_q(buffer2,frame.^2,N);

        if isempty(template)
            windowlen = N;
            %fprintf('N %d buffer %d\n', N, length(buffer));
            [template,startidx] = create_template03(buffer,buffer2,round((p.window+p.windowext*[-1 1])*1e-3*fs),p.rmsthresh,p.template_corr1,p.template_corr2,N,0);
            if 0
             fig = figure(2);clf();
             subplot(211); plot(buffer);
             subplot(212); plot(buffer2);
             drawnow();
            end
            %template = create_template03(buffer,buffer2,round([0.055 0.065]*fs),0.1,0.99,0.9,N);
            if ~isempty(template) && length(template) > (p.window - p.windowext - 1)
                templates{1}(1,:) = template(:)';
                [mxval,mxidx]= max(template);
                %template = buffer(startidx+mxidx-3+[1:numel(template)]);
                template = template(:)';
                tidx = 2;
                spidx = [];
%                nextindex = 1;
                save rn10_template template buffer
            end
        end

        endmark = max(endmark-N,0);

        % Estimate RMS
        RMS = sqrt(mean(buffer.^2));
        %fprintf('RMS:%f\n', RMS);
        if RMS>p.rmsthresh && ~isempty(template)
            % Find best match
            Tsz = size(template,2);
            if isempty(spidx) || (tidx~=spidx)
                %if size(template,1)<tidx | sum(abs(template(tidx,:)))==0,
%                 nzidx = find(sum(abs(template),2)>0);
%                 temp2use = mean(template(setdiff(nzidx,spidx),:),1);
                if size(template,1)< (p.slices+1)
                    temp2use = template(1,:);
                else
                    temp2use = template(p.slices+1,:);
                end
%                 else
%                     temp2use = template(tidx,:);
%                 end
%                 %temp2use = template(1,:);
            else
                temp2use = template(tidx,:);
            end
            c = testmult(buffer,temp2use,N+1,N);

            oldmaxc = maxc2;
            maxc2 = maxc;
            
            oldmidx = midx;

            [maxc,midx] = max(c,[],2);
            if plot_debug
                figure(2);clf();
                plot(buffer(N+1:end));hold on;
                plot(buffer(oldmidx:oldmidx+Tsz-1), 'g');
                plot(0.8 * temp2use, 'r');
                hold off;
                drawnow();
            end
            %fprintf('maxc:%f maxc2:%f cd:%f midx %d\n',maxc, maxc2, (maxc2-maxc + maxc2-oldmaxc), midx);

%             plot(buffer);
%             hold on;plot(temp2use,'r');hold off;
%              cmax_store = [cmax_store,maxc];
             % [tidx oldmaxc maxc2 maxc]
            if maxc<maxc2 && oldmaxc<maxc2  && (maxc2-maxc + maxc2-oldmaxc) > 0.2
                %nextindex = mod(nextindex,p.slices)+1;
                %fprintf('maxc2:%f tidx[%d] nextidx[%d]\n',maxc2,tidx,nextindex);
                idx = oldmidx;
%                 if numel(templates)<tidx | isempty(templates{tidx})
%                     templates{tidx}(1,:) = buffer(idx:idx+Tsz-1)';
%                 else
%                     templates{tidx}(end+1,:) = buffer(idx:idx+Tsz-1)';
%                 end
                if (maxc2>p.corr_remove) | (~isempty(spidx) & (tidx ==spidx)) | ~isempty(spidx)
                    val = buffer(idx:idx+Tsz-1);
                    temp = temp2use';
                    Nfft = 2^nextpow2(Tsz); length(val); %
%                     plot([val,temp]);
                    hf = fft(val-temp,Nfft);
                    ht = fft(temp(:),Nfft);
                    switch p.wt_type,
                        case 0,
                            wt = 0;
                        case 1,
                            wt = 1;
                        case 2,
                            wt = abs(angle(hf.*conj(ht))/pi);
                        case 3,
                            wt = abs([unwrap(angle(hf))-unwrap(angle(ht))]/max(abs([unwrap(angle(hf))-unwrap(angle(ht))])));
                        case 4,
                            wt = 1-abs(linspace(1,0,length(hf)+1)');
                            wt = wt(1:end-1);
                        case 5,
                            wt = min(abs(linspace(1,0,length(hf)+1)'),1-abs(linspace(1,0,length(hf)+1)'));
                            wt = wt(1:end-1);
                        case 6,
                            wt = 1-abs(linspace(-1,1,length(hf)+1)');
                            wt = wt(1:end-1);
                        case 7,
                            wt = 1-abs(linspace(-1,1,length(hf)+1)');
                            wt = wt(1:end-1);
                    end
                    %wt = (wt+fftshift(wt))/2;
                    if p.cutoff > 0,
                        cutoff = round(Nfft*p.cutoff/(fs/2)); % high freq
                        cutoffwt = ones(Nfft,1);
                        cutoffwt(cutoff:(end-cutoff)) = 0;
                        cutoff = round(Nfft*80/(fs/2));        % low freq
                        cutoffwt([1:cutoff,(end-cutoff):end]) = 0;
                    else
                        cutoffwt = 1;
                    end

                    val = real(ifft(max(cutoffwt.*(abs(hf)-p.nc_param*(wt).*abs(ht)),0).*exp(j*angle(hf)),Nfft));
                    val = val(1:Tsz);

                    rmsres2 = sqrt(mean(val.^2));
                    %fprintf('rmsres2: %f\n',rmsres2);
                    rmsres_store = [rmsres_store,rmsres2];
                    if numel(rmsres_store)>30,
                        %p.rmsupdate = nanmean(rmsres_store(1:30))+nanstd(rmsres_store(1:30));
                    end
                    if (rmsres2 < p.rmsupdate)
                        %templates = [templates;buffer(idx:idx+Tsz-1)'];
                        %if size(template,1)<nextindex | (sum(abs(template(nextindex,:)))==0)
                        %    tidx = nextindex;
                        %end
                        if size(template,1)<tidx | (sum(abs(template(tidx,:)))==0)
                            template(tidx,:) = buffer(idx:idx+Tsz-1)';
                            fprintf('adding template rms[%f]tl[%d]fc[%d]\n',rmsres2,tidx,framecount);
                        else
                            template(tidx,:) = p.alpha*template(tidx,:)+(1-p.alpha)*buffer(idx:idx+Tsz-1)';
                            fprintf('updating template rms[%f]tl[%d]fc[%d]\n',rmsres2,tidx,framecount);
                        end
                        if isempty(spidx) | (tidx ~=spidx)
                            template(p.slices+1,:) = p.alpha*temp2use+(1-p.alpha)*buffer(idx:idx+Tsz-1)';
                        end
%                         temp = template(tidx,:)';
%                         hf = fft(val,Nfft);
%                         ht = fft(temp(:),Nfft);
%                         switch p.wt_type,
%                             case 2,
%                                 wt = abs(angle(hf.*conj(ht))/pi);
%                             case 3,
%                                 wt = abs([unwrap(angle(hf))-unwrap(angle(ht))]/max(abs([unwrap(angle(hf))-unwrap(angle(ht))])));
%                             otherwise,
%                         end
                        val(:) = 0;%real(ifft(max(cutoffwt.*(abs(hf)-p.nc_param*(wt).*abs(ht)),0).*exp(j*angle(hf)),Nfft));
                    end
                    buffer(idx:idx+Tsz-1) = val(1:Tsz);
                    %maxc = inf;
                    if idx > endmark,
                        buffer(endmark+1:idx) = 0; %(0.1/RMS)*buffer(endmark+1:idx);
                    end
                    endmark = idx+Tsz-1;
                    if ~isempty(spidx) & (tidx == spidx),
                        [mv,mi] = max(temp2use);
                        buffer(idx+mi) = buffer(idx+mi)+5;
                    end
                elseif isempty(spidx),
                   
                    %tidx = nextindex;
                    if size(template,1)<tidx | (sum(abs(template(tidx,:)))==0)
                        template(tidx,:) = buffer(idx:idx+Tsz-1)';
                        buffer(idx:idx+Tsz-1) = 0;
                        spidx = tidx;
                        %fprintf('adding template tl[%d]fc[%d]\n',tidx,framecount);
                    else
                        template(tidx,:) = p.alpha*template(tidx,:)+(1-p.alpha)*buffer(idx:idx+Tsz-1)';
                        buffer(idx:idx+Tsz-1) = buffer(idx:idx+Tsz-1)-template(tidx,:)';
                        %fprintf('updating template tl[%d]fc[%d]\n',tidx,framecount);
                    end
                end
%                if (size(template,1)==p.slices) | (nextindex == p.slices),
%                    if size(template,1)<p.slices,
%                        template(p.slices,:) = 0;
%                    end
                    tidx = mod(tidx,p.slices)+1;
%                end
            end
        end
    otherwise,
end

function q = push_q(q,sig,N)
q = [q((N+1):end,:);sig];