function [times, Nspk, cumspk ] = synchronydetector(spk, nspk, twin, dt)
%[times Nspk] = synchronydetector(spk, nspk, twin [,dt])


if iscell(spk)   %if several units, concatenate
    spk = mua(spk);
end

if nargin<4
    dt = twin/100;
end
ngauss = round(twin/dt);

[cumspk, allt] = conv_spike(spk, gauss(ngauss,1), dt, 'centered');


overthr = find(cumspk > nspk);

times = [];
Nspk = [];
while ~isempty(overthr)
    i = overthr(1);
    neighb = find(overthr-overthr(1) < twin/dt); %all other times within a window of twin
    [maxo, timo] = max(cumspk(neighb)); %local maximum
    Nspk(end+1) = maxo;
    times(end+1) = allt(i+timo-1);
    
    neighb = find(overthr-overthr(timo) < twin/dt);
    overthr(neighb) = [];   %remove this neighbourhood from the list
    
end

