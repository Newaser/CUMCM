clc;
clear;
global Nodes;
global Actuators;
importData('..\Problems\A\附件1.csv', '..\Problems\A\附件2.csv');
disp(Nodes);
disp(Actuators);

function importData(filepath1, filepath2)
    global Nodes;
    global Actuators;
    
    %import data about Crossed Nodes
    opts = detectImportOptions(filepath1);
    opts.VariableTypes = {'string', 'double', 'double', 'double'};
    opts.DataLines = [2, inf];
    opts.SelectedVariableNames = [1];
    Nodes.ID = readmatrix(filepath1, opts);
    opts.SelectedVariableNames = [2:4]; 
    Nodes.Pos = readmatrix(filepath1, opts);

    %import data about Actuators
    opts = detectImportOptions(filepath2);
    opts.VariableTypes = {'string', 'double', 'double', 'double', ...
        'double', 'double', 'double'};
    opts.DataLines = [2, inf];
    opts.SelectedVariableNames = [1];
    Actuators.ID = readmatrix(filepath2, opts);
    opts.SelectedVariableNames = [2:4];
    Actuators.BottomPos = readmatrix(filepath2, opts);
    opts.SelectedVariableNames = [5:7];
    Actuators.TopPos = readmatrix(filepath2, opts);

    % disp(Nodes);
    % disp(Actuators);
    % disp(Nodes.ID);
    % disp(Nodes.Pos);
    % disp(Nodes.ID);
    % disp(Actuators.ID);
    % disp(Actuators.BottomPos);
    % disp(Actuators.TopPos);
end
