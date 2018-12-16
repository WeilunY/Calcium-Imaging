function [ rerun ] = MC_gui_final( )
%INITIALGUI Summary of this function goes here
%   Detailed explanation goes here
% Help goes here.

% Code running parameters
rerun = false;

guiHeight = 200;
guiWidth = 300;

%%% POSITION: [from-left from-bottom width height]

% 'position', all distance units in pixels: [xdist ydist width height]
S.fh = figure('units','pixels',...
            'position',[300 300 guiWidth guiHeight],...
         	'menubar','none',...
          	'name','Ca Code',...
         	'numbertitle','off',...
        	'resize','off');
S.title = uicontrol('Style','text','fontsize',16,...
            'fontweight','bold',...
            'String','Rerun Motion Correction?',...
            'Position',[5 guiHeight-60 guiWidth 50]);
S.horozontalLine = uicontrol('Style','text','fontsize',12,...
            'fontweight','bold',...
            'String','',...
            'Position',[0 guiHeight-55 600 2],'backgroundcolor','k');
%%% START Button
S.buttonRerun = uicontrol('style','push',...
            'unit','pix',...
         	'position',[10 75 280 50],...
            'fontsize',12,...
            'fontweight','bold',... 
          	'string','RERUN MC',...
           	'callback',@button_rerun_call);
S.buttonContinue = uicontrol('style','push',...
            'unit','pix',...
         	'position',[10 15 280 50],...
            'fontsize',12,...
            'fontweight','bold',... 
          	'string','CONTINUE',...
           	'callback',@button_continue_call);


guidata(S.fh,S)
movegui('center')

function [] = button_continue_call(varargin)
    % Callback for pushbutton.
    S = guidata(gcbf);
    
    rerun = false;
    
    close all force
end

function [] = button_rerun_call(varargin)
    % Callback for pushbutton.
    S = guidata(gcbf);
    
    rerun = true;
    
    close all force
end

function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
%   NOT USED, AL DONE IN ABOVE FUNCTION. LEFT FOR REFERENCE
end

%%% Pauses the program until the figure 'S.fh' closes
waitfor(S.fh)

end

