function [  ] = dummyDataGenerator( )
%DUMMYDATAGENERATOR dummyDataGenerator(1)
%   Detailed explanation goes here
frames = 1143; % exactly 5 minutes of frames
neurons = 10;
data = ones(neurons, frames);

x = linspace(0,300,300)
plot(x, 0.00004*x.^2)

pause

rng default  %initialize random number generator
rng('shuffle')
for cell=1:neurons
    %noise1 = 0.3*rand(1,frames);
    %noise2 = 0.15*rand(1,frames);
    noise = 0.15*randn(1,frames); % White Gaussian Noise > flat distribution
    actionPots = zeros(1,frames);
    for i=1:5
        index = round(frames*rand(1,1));
        modelExp = getRandExp();
        actionPots(index:index+3) = modelExp;
    end
   % data(cell,:) = actionPots; % no noise
   % data(cell,:) = noise1 + noise2 + actionPots; % flat distribution
    data(cell,:) = noise + actionPots;
end

hold off
inc = 1;
for cell=1:neurons
    plot( data(cell,:)+(inc-1),'linewidth',1.1)
    inc = inc+1;
    hold on;
end

assignin('base', 'fakeData', data)

end
% 2 big 2 mid 2 s
function [modelExp] = getRandExp()
% Returns a vector 4 units long of a pseudo exponential
    expMatrix = [[1, 0.4, 0.25, 0.19]
                 [1, 0.7, 0.6, 0.4]
                 [1, 0.24, 0.07, 0]
                 [1, 0.85, 0.75, 0.65]];
    size(expMatrix);
    rng default 
    rng('shuffle')
    indexRand = ceil(4*rand(1,1));
    modelExp = expMatrix( indexRand,: );
    % std of noise to some known level of this scale %%%%%%
    eventSizes = [0.8, 0.55, 0.3];
    eventSizes = [1 .8 .7];
    amplitude = eventSizes(ceil(3*rand(1,1)));
    modelExp = modelExp*amplitude;
end