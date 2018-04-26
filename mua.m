function spk = mua(spikes)
%spk = mua(spikes)
%spikes is a cell array (each cell for a single neuron) of spike timing arrays
%spk is a spike timing array concatenating all spike times form the
%different spikes

spk = [];
for i=1:length(spikes)
    spk = [spk spikes{i}];
end
spk = sort(spk);