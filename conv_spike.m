function [y times] = conv_spike(spk, kern, dt, lag, trange)
%[y times] = conv_spike(spk, kern, dt [,delag], [trange])
%
%-spk : train of impulses/spikes (each value is the timing of a spike)
%-kern : the kernel by which it must be convoluted
%-dt : time resolution  of the kernel
%-lag : onset of the kernel (must be negative to avoid pahse-shifting)
%      use c'entered' to center the kernel on 0
%-trange : time range of output vector

%warning ; there might be a small time-shift (dt/2)

%see also smooth_spike (not written yet hehe)



if nargin < 4,
    lag = dt;
elseif same(lag, 'centered'),
    lag =  - length(kern)*dt/2;
end

%if cell array : repeat function for each element
if iscell(spk),
    if nargin<5
    trange = [];
    end
    vararg = {kern, dt, lag, trange};
    vararg = vararg(1:nargin-1);
    y = cell(size(spk));
    for i = 1:numel(spk)
        [y{i} times] = conv_spike(spk{i}, vararg{:});
    end
    return;
    
end


if nargin <5,
    
    if isempty(spk),
        y = [];
        times = [];
    end
    minspk = min(spk);
    spk = spk - minspk;
    maxspk = ceil(max(spk/dt));
    shape = 'full';
    times = minspk + dt*[0:maxspk+length(kern)-1] + lag;
    
else   %specified time range
    minspk = trange(1);
    maxspk = ceil((trange(2)-trange(1))/dt);
    shape = 'same';
    times = minspk + dt*[0:maxspk] + lag + dt*length(kern)/2;
    
    if isempty(spk),
        y = zeros(1,length(times));
        return;
    end
    spk = spk - trange(1);
    
    
end

ispk = 1+round(spk/dt);
ispk( ispk<0 ) = [];   %if using trange, remove spikes outside of the range
ispk( ispk>maxspk+1) = [];

ccc = zeros(1,maxspk+1);
for i = 1:length(ispk),
    ccc(ispk(i)) = ccc(ispk(i))+1;
end

y = conv(ccc, kern, shape);






