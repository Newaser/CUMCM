clc;
clear;
global Nodes;
global Actuators;
global Reflectors;
importData('..\Problems\A\附件1.csv', '..\Problems\A\附件2.csv', ...
    '..\Problems\A\附件3.csv');
disp(Nodes);
disp(Actuators);
disp(Reflectors);

function importData(filepath1, filepath2, filepath3)
    global Nodes;
    global Actuators;
    global Reflectors;
    
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
    
    %import data about Reflector elements
    opts = detectImportOptions(filepath3);
    opts.VariableTypes = {'string', 'string', 'string'};
    opts.DataLines = [2, inf];
    opts.SelectedVariableNames = [1:3];
    Reflectors = readmatrix(filepath3, opts);

    % disp(Nodes);
    % disp(Actuators);
    % disp(Reflectors);
    % disp(Nodes.ID);
    % disp(Nodes.Pos);
    % disp(Nodes.ID);
    % disp(Actuators.ID);
    % disp(Actuators.BottomPos);
    % disp(Actuators.TopPos);
end
