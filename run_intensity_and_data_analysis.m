function [  ] = run_intensity_and_data_analysis( fileData )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    load([fileData,' motionCut_Output.mat']);
    load([fileData,' segOutput.mat']);
    analyzeIntensityGUI( intavgNPCorr, [fileData,' intensity data']);
    load([fileData,' intensity data.mat']);
    intensityZ( stackAdjustedCut, intavgNPCorr, firingTimes, binaryFiring, firedNeurons, fileData);

end

