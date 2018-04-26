function rasterplot(spikes, block, style)
%rasterplot(spikes) plots a raster plot where spikes is a cell where the
%n-th element represents the spiking times of the n-th neuron
%
%rasterplot(spikes, block)

N = length(spikes); %number of neurons

if nargin ==1,
    block = N;
    style = {'k'};
elseif nargin ==2,
    style = { 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'};
end

ish = ishold;

sblock = cumsum(block);
if sblock(end) ~= N,
    error('the sum of neurons in each block should equal to the number of neurons'); 
end

for i=1:N,
    if ~isempty(spikes{i})
        ss = find(sblock>=i, 1);
        if isnumeric(style{ss})
                    plot([spikes{i} ; spikes{i}], i-1:i, 'Color', style{ss});

        else
        plot([spikes{i} ; spikes{i}], i-1:i, style{ss});
        end
    end
    if i==2,
        hold on;
    end
end

if ~ish,
    hold off;
end

