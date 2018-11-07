function [  ] = run_segmentation( fileData,region )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    load([fileData,' motionCut_Output.mat']);
    tif2IntensitiesMatt(stackAdjustedCut,10,[fileData,' segOutput'], region);

    if SG_gui_final
        run_segmentation( fileData,region )
    end
end

