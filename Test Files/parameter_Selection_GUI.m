function [ brainRegion ] = parameter_Selection_GUI(  )
%PARAMETER_SELECTION_GUI Allows user to alter functionality
%   Changes certain variables to suit the user

brainRegion = 'Test';

% This is the figure we will put everything onto
f = figure('Visible','off');
% Gives the figure a name and removes annoying numbering
figure('Name','Selection','NumberTitle','off');

pos = get(gcf, 'Position'); %// gives (x left, y bottom, width, height)
width = pos(3);
height = pos(4);

% Position, [right up width height]     <- format of the Position option
% h stands for handle, allows us to check the values of each checkbox
h1 = uicontrol('Style','checkbox','Callback',@checkbox1_Callback,...
  'Position',[ 50 height-50 50 50],'String', 'CA1','Value',1 );
h2 = uicontrol('Style','checkbox','Callback',@checkbox2_Callback,...
  'Position',[ 150 height-50 50 50],'String', 'DG','Value',0 ); 
% get(h1, 'Value'); returns the value of whatever is associated with h1
% h1.Value = 0;  simply changes the value of what is associated with h1

% This function is called when checkbox1 is toggled
function checkbox1_Callback(hObject, eventdata, handles)
    % get(hObject,'Value')  returns toggle state of checkbox1
    if (get(hObject,'Value') == get(hObject,'Max'))
        h2.Value = 0;
        brainRegion = 'CA1';
    else
        h2.Value = 1;
    end
end
% This function is called when checkbox2 is toggled
function checkbox2_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value') == get(hObject,'Max'))
        h1.Value = 0;
        brainRegion = 'DG';
    else
        h1.Value = 1;
    end
end
% pause(); will stop the code here until the user types a character below
%   indicating that they have finished selecting
pause();
% Prints the value set by the user
disp(brainRegion);

end