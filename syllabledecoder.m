function [putativeonsets,  yy, LFPt] = syllabledecoder(filename, datatype, paramfile, doplot, domex)
% putativeonsets = syllabledecoder(wavfile, datatype)
% uses theta neural network to decode syllable boundaries from WAV
%putativeonsets is given is ms
%
%putativeonsets = syllabledecoder(filename, datatype, paramfile)
%paramfile is a .MAT with parameters for theta neural network, spectral filter and temporal convoluution filter applied on auditory channels
%
%
%putativeonsets = syllabledecoder(filename, datatype, paramfile,  1)
%to plot auditory channels, spikes and putative onsets
%
%[putativeonsets  Vmean] = syllabledecoder(...)
%Vmean is the average membrane potential of excitatory neurons - continuous
%variables sampled at dt = 0.1 ms(continuous except for reset after
%spikes)
%
%[putativeonsets, Vmean, LFP] = syllabledecoder(...)
% LFP trace
%
% Example:
%[putativeonsets, Vmean, LFP] = syllabledecoder('SI711.wav','native','speechparser_param',1);

if nargin <4
    doplot = 0;
end
if nargin<5
    domex = 0; %control whether uses the MEX file for network simulation
end

inptype = 'channel';  %'envelope', 'channel'
fe = 1000;        %sampling frequency of downsampled channel/envelope (Hz)

%% open WAV file
[x, wav_fe]=audioread(filename,datatype);
% fid2 = fopen(filename, 'r');
% fseek(fid2,1024,-1); % Skip past header, which is 1024 bytes long, for TIMIT sentences
% x=fread(fid2,inf,'int16'); % 16-bit signed integers
% x = x(:)'; % Row vector
% fclose(fid2);

if strcmp(inptype, 'channel')
    %% process WAV through Chi and Shamma's cochlear model
    
    %nb=1;
    octshft=-1; % -1 for 8kHz
    rf = 8000;    % resample at 8 kHz
    %DBnoise = 100;
    
    loadload;  % Loads colormap, filterbank and parameters
    %warning('off', 'stats:categorical:subsasgn:NewLevelsAdded');
    
    % downsample
    yr = resample(double(x), rf, wav_fe);
    
    % time, frequency
    %L = length(yr);
    paras(4) = octshft;							% octave shift
    paras(1:2) = paras(1:2) / 2 ^(octshft+1);	% frame size and step
    paras(3) = -2;								% linear
    %SF = 16000 * 2^octshft;						% sample rate
    
    % waveform
    yr = unitseq(yr);                 % sequence normalization: on centre et on rÈduit
    F = 1/2*(440 * 2 .^((-31:97)/24));  %center frequencies of all channels
    
    % tranform wav into auditory spectrogram
    AS = wav2aud(yr, paras);
    if 0 % change to 1 to plot auditory spectrogram
        subplot(2,1,2);
        aud_plot(AS, paras); % plot auditory spectrogram (log frequency scale)
    end
    
    %upsample
    l = size(AS,1);
    x = AS(  ceil(fe/1000/paras(1):fe/1000/paras(1):l)  ,:);
    
else
    %% extract the envelope
    
    %warning('off', 'stats:categorical:subsasgn:NewLevelsAdded');
    
    %envelope is just absolute value of hilbert transform
    h = hilbert(x);
    x = abs(h);
    
    %downsample
    x = smooth(x, wav_fe/fe);
    x = downsample(x, wav_fe/fe);
    
end


%% PARAMETERS FOR THETA NEURAL NETWORK

Tstart = 500; %onset of sentence stim (after 500 ms)

dt = .1;     % time step for Euler equation
dtt = 5;       %temporal resolution of spectrograms (ms)

%colorz = {[68 120 0]/256,[160 202 74]/256,[0 119 162]/256, [57 183 216]/256}; %colors for D, S, F and E neurons resp.


nD = 10;    %number of slow excitatory neurons
nS = 10;    % number of slow inhibitory neurons
nchan = 32; %auditory channels

%dem = load(chanfilterfile);
%Ichannel2D = dem.sta_chan;

load(paramfile);


n =  nD + nS;
DD = 1:nD;
SS = n+(1-nS:0);
II = SS;

%%%%%%%%%%%%%%%%%%
C = 1;   %membrane capacitance density (F/cm2)

VL = -67;
VI = -80;
VE = 0;
Vsyn = zeros(n);
Vsyn(DD,:) = VE;
Vsyn(II,:) = VI; %[VE*ones(nE+nD,n) ; VI*ones(nI,n)];

gL = .1;%conductance density (mS/cm2)

C = C*ones(1,n);
%VLs = -70;
Vthr = -40*ones(1,n);
Vres = -87*ones(1,n);
VL = VL*ones(1,n);

gLs = 0.02;
gL = [gLs*ones(1,nD) gL*ones(1,nS)];
gL(DD) = gLD*gL(DD);

%%%%%%%%%%%%%%%% CONNECTIVITY

g = zeros(n);
g(DD,DD) =  0/nD;
g(SS,SS) = 5/nS ;
g(DD,SS) = 50/nD;
g(SS,DD) = 1.5/nS;
frfr = 20;
g(DD,:) = g(DD,:)/frfr;
g(DD,SS) = gDS * g(DD,SS);
g(SS,DD) = gSD * g(SS,DD);
g(SS,SS) = gSS * g(SS,SS);

%synpatic time constants
tauR = zeros(n,1);
tauD = zeros(n,1);

tauR(DD) = .2;  %synaptic rise time constant
tauR(SS) = 5;
tauR(DD) = frfr*tauR(DD);

tauD(DD) = 50; % synaptic decay time constant
tauD(SS) = 25;
tauD(DD) = tauDD * tauD(DD);
tauD(SS) = tauDS * tauD(SS);


%% %%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
%membrane potentials (in mV)
%load membrane potential distribution for gamma neurons (randomly distrib
%around the spiking cycle

V = VL + (rand(1,n)-.5).*(Vthr-VL);
V(DD) = -80+40*rand(1,nD);

ss = zeros(1,n);


%% %%%%%%%%%%%%%%% INPUTS %%%%%%%%%

%current density (in ï¿½A/cm2)
Idc = zeros(1,n);
Idc(SS) = 0.1;
Idc(DD) = .6;
Idc(SS) = IdcS * Idc(SS);
Idc(DD) = IdcD * Idc(DD);

%Inoise
knoise = zeros(1,n);
knoise(DD) = 2*knoizD*knoiz;
knoise(SS) = 2*knoiz;

Iext = zeros(1,n);
Iext(DD) = 1.5;
Iext(SS) = 0;

Tons = Tstart;
%Toff = Tons + stidur';
%Tinp = [Tons Toff];
%Tinpstep = ceil(Tinp/dt);
%It = zeros(1,TT);

switch inptype
    case 'envelope'
        envnorm = 3000; %input normalization factor
        
        %load data from envelope file and syllable time file
        
        
        env = x/envnorm; %normalize
        stidur = ceil(length(env)/(fe/1000)/dt)*dt; %duration of sentence (ms)
        ups = dt*(fe/1000);                %upsampling
        env = env(ceil(ups:ups:stidur*(fe/1000)));%resample at dt
        
        T = Tstart + stidur + 200; % Tstart ms before sentence, 200 ms after
        T = ceil(T/dtt)*dtt;
        TT = ceil(T/dt);
        It = zeros(1,TT);
        It(ceil(Tstart/dt)+(1:length(env))) = env;  %and start at Tstart
        
        Tons = Tstart;
        Toff = Tons + stidur;
        Tinp = [Tons Toff];
        
    case 'channel'
        
        %normalisation
        chan = x'/3; %maximum over sentences is around 60
        whichchan = (1:nchan)*floor(128/nchan); %each E cell receives input from a different channel from cochlea
        stidur = ceil(size(chan,2)/(fe/1000)/dt)*dt;
        ups = dt*(fe/1000);
        
        T = Tstart + max( stidur+200, 2000); % Tstart ms before sentence, 200 ms after (at least 2000)
        T = ceil(T/dtt)*dtt;
        TT = ceil(T/dt);
        
        IIt = zeros(n,TT);
        TTchan = stidur*(fe/1000);
        
        %spectral filter
        chaninput = Ichannel2D*chan(whichchan, 1:TTchan);
        
        % temporal convolution
        %         Dconv = load(tempfilterfile);
        %        if ~isempty(Dconv)
        %             if ischar(Dconv)
        %                 Dconv = load(Dconv);
        %             end
        %            chaninput = downsample( smooth(chaninput,fe/Dconv.fe), fe/Dconv.fe)';
        %            chaninput = filter(Dconv.Dconv, [1 zeros(1,length(Dconv.Dconv)-1)], chaninput);
        %            chaninput = chaninput(ceil(Dconv.fe/fe:Dconv.fe/fe:length(chaninput)));
        %        end
        chaninput = downsample( smooth(chaninput,fe/conv_fs), fe/conv_fs)'; % smooth and downsample signal to sampling rate of temporal kernel (100 Hz)
        chaninput = filter(Dconv, [1 zeros(1,length(Dconv)-1)], chaninput); % convolve with temporal kernel
        chaninput = chaninput(ceil(conv_fs/fe:conv_fs/fe:length(chaninput))); % upsample back to 1000 Hz
        IIt(DD, ceil(Tstart/dt)+(1:stidur/ups)) = Iext(DD)' * chaninput(:,ceil(ups:ups:TTchan));
        
        Tons = Tstart;
        Toff = Tons + stidur;
        % ISI = [];
        % nbinp = 1;
        % Tinp = [Tons Toff];
        %   warp = 1;size(IIt)
        
        
end

IIt(:,TT+1:end) = [];
It = IIt(DD(1),:)/Iext(DD(1));


%% %%%%%%%%%%%%%%%%%%%%% NEURAL NETWORK SIMULATION %%%%%%%%%%%%%%%

if domex  %simulate with MEX file
    
    Vsyn = Vsyn(:,1);
    
    %   if 1
    [spikes, Vt, LFPt] = lif2(n, dt, TT, V, ss', ...
        C, VL, Vthr, Vres, Vsyn', ...
        gL, g, tauR', ...
        tauD', Idc, IIt + repmat(knoise'/sqrt(dt),1,TT).*randn(n,TT));
    %   else
    %       [spikes ] = lif3( TT); % just for debugging -> remove that
    %   end
    %    nspike = cellfun(@length, spikes);
    %    disp([length(spikes) sum(nspike)]);
    
else
    
    %   IL = zeros(1, n);
    %   Io = zeros(1, n);
    %   Isyn = zeros(1,n);
    %   Im = zeros(1,n);
    XXXt = zeros(TT,n);
    
    %  winf = zeros(1,n);
    %  taum = zeros(1,n);
    
    spikes = cell(1,n);
    nspike = zeros(1,n);
    %  suprathr = zeros(1,n);
    last = zeros(1,n);
    % nuspk = zeros(1,n);
    
    ss = zeros(n,1); %not n x n, because the activation is equal for all synapses coming from the same neuron
    sr = zeros(n,1); %rise component
    
    Vt = zeros(TT, n); %membrane potentials history
    LFPt = zeros(1,TT);
    
    
    for t = 1:TT
        
        %standard leak current
        IL = gL .*(VL-V);
        
        %M-current
        %         winf = 1./(1+exp(-(V+35)/10));
        %         %taum = 400./ (3.3*exp((V+35)/20) + exp(-(V+35)/20));
        %         taum = 100./ (3.3*exp((V+35)/20) + exp(-(V+35)/20));
        %         w = w + dt*(winf-w)./taum;
        %         Im = gM .* w .*(VK - V);
        Im = zeros(1,n);
        
        %external currents
        Io = Idc + It(t)*Iext;
        
        %synaptic current
        %sr = 1/2 * (1 + tanh(V'/10)) .* (1-ss) ./tauR;
        sr = sr .* (1- dt./tauR);
        ss = ss + dt* (sr - ss)./tauD;
        Iss = g .* repmat(ss,1,n) .* (Vsyn-repmat(V,n,1));
        Isyn = sum(Iss);
        
        %internal noise
        Inoise = knoise .* randn(1,n) / sqrt(dt);
        
        %membrane potential
        V = V + dt * (IL + Im + Io + Isyn + Inoise)./C;
        
        %LFP
        LFP = - sum(sum(abs(Iss(:,DD))));
        
        if any(isnan(V))
            error('membrane potential reaches nan value');
        end
        
        %V(nws+1:n) = VL;
        
        suprathr = (V>=Vthr) ;
        for i= find(suprathr)
            nspike(i) = nspike(i)+1;
            spikes{i}(nspike(i)) = t*dt;
            last(i) = t*dt;
            V(i) = Vres(i);
            sr(i) = sr(i)+1;            
        end
        last = suprathr;
        XXXt(t,:) = ss;
        Vt(t,:) = V;
        LFPt(t) = LFP;
    end
    
end

%% %%%%%%%%%%%%%%%%  ANALYSIS %%%%%%%%%%%%%%%%%%%

% sum of membrane potentials over excitatory neurons, starting at wav
% onset
yy = sum(Vt(DD,ceil(Tstart/dt)+1:end),1);



%detect neuron synchrony
[synS, nsynS]= synchronydetector(spikes(SS), nS/5,15,.03);
[synD, nsynD]= synchronydetector(spikes(DD), nD/5,15,.03);
%periodS = diff(synS);


%% syllable analysis
putativeonsets = synS - Tons;

%% %%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%

if doplot %||1
    
    %% raster plot
    figure(); clf;
    subplot('Position', [0 .85 1 .15]); hold on;   % LFP
    plot((dt:dt:T)-Tstart, LFPt, 'k');
    %plot((dtt:dtt:T)-Tstart, thetaband, 'Color', colorz{1}, 'Linewidth', 2);
    
    xlim([-1-Tstart T+1-Tstart]); axis off
    
    plotchan = strcmp(inptype, 'channel');
    
    subplot('Position', [0 0+.4*plotchan 1 .8-.4*plotchan]);   %raster plot
    hold on;
    spkshift = cellfun(@(x) x-Tstart,spikes([DD SS  ]),'Un',false); % shift spike times
    rasterplot(spkshift,[nD nS  ],  {'r','k'});hold on; % plot raster for Te and Ti neurons
    
    %stim
    if strcmp(inptype, 'envelope')
        rasterplot({sylon syloff}, [], 'b')
    end
    axis([-1-Tstart T+1-Tstart 0 n])
    if ~plotchan
        plot((dtt:dtt:T)-Tstart, 10*smooth(Itt,5/dtt),'r')
    else
        axis off;
        subplot('Position', [0 0 1 .35]);
        imagesc([0 Toff-Tons], [1 size(chan,1)], chan(:,ceil(ups:ups:stidur*fe/1000)));  %plot aud spectrogram
        Nchanplot = size(chan,1);
        
        colormap(flipud(gray));
        axis([-1-Tstart T+1-Tstart 0 Nchanplot]);
        axis xy; axis off;
        %         if strcmp(inptype, 'channel')
        %             hold on;
        %             plot(repmat(sylon,2,1),[0 Nchanplot], 'Color', colorz{1});
        %         end
    end
end


end