clc;
clear;
global Nodes;
global Actuators;
global Reflectors;
importData('..\Problems\A\附件1.csv', '..\Problems\A\附件2.csv', ...
    '..\Problems\A\附件3.csv');
% disp(Nodes);
% disp(Actuators);
% disp(Reflectors);

%pre-definition
node_num = length(Nodes.ID);
R = 300;
F = 0.466*R;
r = R-F;
alpha = 0;
beta = 90;
r_cabin = 0.5;

hold on
%draw all the nodes
plot3(Nodes.Pos(1:end, 1), Nodes.Pos(1:end, 2), Nodes.Pos(1:end, 3), ...
    '.k', 'markersize', 8);

%draw the center point
plot3(0, 0, 0, 'or', 'markersize', 2);

%draw the light source
[s1, s2, s3] = drawSource(alpha, beta, 200, 1.5);
S = [s1, s2, s3];

%draw all the Actuators
for i = 1:node_num
    x = [Actuators.TopPos(i, 1), Actuators.BottomPos(i, 1)];
    y = [Actuators.TopPos(i, 2), Actuators.BottomPos(i, 2)];
    z = [Actuators.TopPos(i, 3), Actuators.BottomPos(i, 3)];
    line(x, y, z, 'color', '#333333', 'LineWidth', 5);
end

%draw all the Tie-down cables
for i = 1:node_num
    x = [Actuators.TopPos(i, 1), Nodes.Pos(i, 1)];
    y = [Actuators.TopPos(i, 2), Nodes.Pos(i, 2)];
    z = [Actuators.TopPos(i, 3), Nodes.Pos(i, 3)];
    line(x, y, z, 'color', '#A2A2A2');
end

%draw all the reflectors
for i = 1:length(Reflectors)
    vertexes = zeros(3,3);
    for j = 1:node_num
        for k = 1:3
            if strcmp(Reflectors(i, k), Nodes.ID(j))
                vertexes(k,:) = Nodes.Pos(j, :);
            end 
        end
    end
    drawTriangle(vertexes);
end

%draw focal hemisphere
% drawFocalHemi(r);

%draw feedback cabin
drawFeedbackCabin(S, r, r_cabin);

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

function drawTriangle(vertexes)
    for i = 1:3
        x = [vertexes(i,1), vertexes(rem(i, 3) + 1,1)];
        y = [vertexes(i,2), vertexes(rem(i, 3) + 1,2)];
        z = [vertexes(i,3), vertexes(rem(i, 3) + 1,3)];
        line(x, y, z, 'color', 'k');
    end
end

function [x, y, z] = drawSource(a, b, r, k)
    a_ = a*pi/180;
    b_ = b*pi/180;
    [x, y, z] = sph2cart(a_, b_, r);
    plot3(x, y, z, '.y', 'markersize', 15);
    line([x -k*x], [y -k*y], [z, -k*z],'color', 'y', 'LineStyle', '--');
end

function drawFocalHemi(r)
    [x, y, z] = sphere;
    x = x*r;
    y = y*r;
    z = z*r;
    meshc(x,y,z);
end

function drawFeedbackCabin(S, r, r_cabin)
    syms x y z k
    assume(x<=0 & y<=0 & z<=0);
    eqn1 = x == k*S(1);
    eqn2 = y == k*S(2);
    eqn3 = z == k*S(3);
    eqn4 = x^2 + y^2 + z^2 == r^2;
    Sol = solve([eqn1 eqn2 eqn3 eqn4], [x y z k], 'ReturnConditions', true);
    drawCircle(S, r_cabin, [Sol.x, Sol.y, Sol.z]);
end

function drawCircle(n, r, c)
    theta=(0:2*pi/100:2*pi)'; %theta角从0到2*pi
    a=cross(n,[1 0 0]); %n与i叉乘，求取a向量
    if ~any(a) %如果a为零向量，将n与j叉乘
        a=cross(n,[0 1 0]);
    end
    b=cross(n,a); %求取b向量
    a=a/norm(a); %单位化a向量
    b=b/norm(b); %单位化b向量

    c1=c(1)*ones(size(theta,1),1);
    c2=c(2)*ones(size(theta,1),1);
    c3=c(3)*ones(size(theta,1),1);

    x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%圆上各点的x坐标
    y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%圆上各点的y坐标
    z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%圆上各点的z坐标

    plot3(x, y, z,'color', 'b', 'LineWidth', 3);
end

