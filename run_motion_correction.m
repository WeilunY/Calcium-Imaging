function [  ] = run_motion_correction( fileData,PathName,FileName1,FileName2,motionData )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [performMotionCorrection, line_to_line_motion_correction] = MC_gui_initial();
    
    globalMC(FileName1,PathName,FileName2,PathName,motionData, performMotionCorrection, line_to_line_motion_correction);
    load([fileData,'.mat']);
    if motionData
        load([PathName,FileName2]);
    else
        Y1 = 0;
    end
        %%%%% Run Motion Removal To Clean Up Movie
    motionCutterGUI(motionCompensation, stackAdjusted,Y1,[fileData,' motionCut_Output'],motionData);
    load([fileData,' motionCut_Output.mat']);
        %%%%% Save The Final Tiff Stack For Reference Later
    toTiff(stackAdjustedCut,[fileData,' CUT']);
    
    rerun = MC_gui_final();
    if rerun
        run_motion_correction()
    end

end

