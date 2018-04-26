function spk = spike_shift(spk, T, match)
%SPK = spike_shift(SPK, T) shift spike times in spike data by time T.
%SPK can be in any supported format for spike data.
%
%SPK = spike_shift(SPK, T, 'neuron') allows to use one shift for each neuron in spike data SPK.
% T must be a vector of length n_neuron, where element n corresponds to shift for neuron n.
%
%SPK = spike_shift(SPK, T, 'trial') allows to use one shift for each neuron in spike data SPK.
% T must be a vector of length n_neuron, where element n corresponds to shift for neuron n.
%
% See also SPIKE_WINDOW, SPIKE_FORMAT, SPIKE_SPIKETIMELOCK

if nargin<3,
    match = 'none';
elseif ~any(strcmpi(match, {'trial','neuron'})),
    error('third argument must be either ''trial'' or ''neuron''');
end

form = spike_format(spk);

% number of neurons and trials in spike data
[n_neuron,n_trial] = spike_properties(spk);

% make sure T has the right size
switch match,
    case 'none',
        if length(T)>1,
            error('T must be a scalar');
        end
    case 'neuron'
        if length(T)~=n_neuron,
            error('with option ''trial'', T must be a vector of length n_neuron');
        end
    case 'trial'
        if length(T)~=n_trial,
            error('with option ''trial'', T must be a vector of length n_trial');
        end
end

% shift values
switch form
    
    % cell format
    case 'cell',
        % build simply a n_trial x n_matrix for T
        switch match,
            case 'none',
                T = T*ones(n_trial, n_neuron);
            case 'neuron',
                T = repmat(T(:)',n_trial,1);
            case 'trial',
                T = repmat(T(:),1,n_neuron);
        end
        for n =1:n_neuron
            for t=1:n_trial
                spk{t,n} = spk{t,n} + T(t,n);
            end
        end
        
        % numeric or struc formats
    case {'numeric','struct'},
        spk = spike_format(spk, 'struct');
        switch match,
            case 'none',
                spk.time = spk.time + T;
            case 'neuron',
                for n=1:n_neuron,
                    thisn = spk.neuron ==n;
                    spk.time(thisn) = spk.time(thisn) + T(n);
                end
            case 'trial',
                for t=1:n_trial,
                    thist = spk.trial ==t;
                    spk.time(thist) = spk.time(thist) + T(t);
                end
        end
        
        spk = spike_format(spk, form); % format back to numeric if this was the case
        
    case 'incorrect',
        error('SPK format is incorrect. Type "help spike_format" for more information on spike formats');
end