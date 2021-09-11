clc; clear;
%% DATA IMPORT PART
% Import the components
global Nodes;
global Actuators;
global Reflectors;
importData('..\Problems\A\附件1.csv', '..\Problems\A\附件2.csv', ...
    '..\Problems\A\附件3.csv');

%% PRETREATMENT PART
% Pre-definition
R = 300;
R_FAST = 500*0.5;
F = 0.466*R;
r = R-F;
% alpha_degree = 0;
% beta_degree = 90;
alpha_degree = 36.795;
beta_degree = 78.169;
r_cabin = 1*0.5;
sourceDist = 200;

% Pre-caculate
node_num = length(Nodes.ID);
alpha = alpha_degree*pi/180;
beta = beta_degree*pi/180;

%% Coordinate System Transformation Part
[alpha, beta] = CoordSysTrans(alpha, beta, node_num);
    

%% Caculation Part
% About FAST
    % getFASTSphCenter();
    % getFASTCaliberCenter();
    % getFASTCaliberCircle();
    [X_VS, Y_VS, Z_VS] = getVirtualSphere(R, R_FAST);

% About Light Source
    [x_S, y_S, z_S] = getSourcePoint(alpha, beta, sourceDist);
    [x_Proj, y_Proj, z_Proj] = getProjectionPoint([x_S, y_S, z_S], R);

% About Feedback Cabin
    [x_cabC, y_cabC, z_cabC] = getCabinCenter([x_S, y_S, z_S], r);
    [X_cabCir, Y_cabCir, Z_cabCir] = ...
        getCabinCircle([x_cabC, y_cabC, z_cabC], [x_S, y_S, z_S], r_cabin);
    [X_cabD, Y_cabD, Z_cabD] = ...
        getCabinDisk([x_cabC, y_cabC, z_cabC], alpha, beta, r_cabin);

% About Paraboloid
    [x_PClbC, y_PClbC, z_PClbC] = getParaCaliberCenter(alpha, beta, R);
    [X_PClbCir, Y_PClbCir, Z_PClbCir] = ...
        getParaCaliberCircle([x_PClbC, y_PClbC, z_PClbC], R/2);
    [X_optP, Y_optP, Z_optP] = ...
        getOptPara(R, F, [x_PClbC, y_PClbC, z_PClbC]);


%% Graphic Plot Part
hold on
% About Components
    drawNodes();
    drawActuators(node_num);
    drawTiedownCables(node_num);
    % drawReflectors(node_num); % Warning: ...

% About FAST
    drawFASTSphCenter();
    % drawFASTCaliberCenter();
    drawFASTCaliberCircle();
    % drawVirtualSphere(X_VS, Y_VS, Z_VS);

% About Light Source
    drawSourcePoint(x_S, y_S, z_S);
    % drawProjectionPoint();
    drawLightPath([x_S, y_S, z_S], [x_Proj, y_Proj, z_Proj]);

% About Feedback Cabin
    % drawCabinCenter();
    drawCabinCircle(X_cabCir, Y_cabCir, Z_cabCir);
    % drawCabinDisk(X_cabD, Y_cabD, Z_cabD);    % Bad Method
    
% About Paraboloid
    % drawParaCaliberCenter();
    drawParaCaliberCircle(X_PClbCir, Y_PClbCir, Z_PClbCir);
    drawOptPara(X_optP, Y_optP, Z_optP);
    
%% FUNCTION PART
%% About Coordinate System
function [a, b] = CoordSysTrans(a, b, node_num)
    global Nodes;
    global Actuators;
    
    % Transformation Matrix A:
    A = [ cos(a),	sin(a),      0;
         -sin(a),	cos(a),      0;
               0,	     0,      1];
    % Transformation Matrix B:
    B = [ sin(b),        0, -cos(b);
               0,	     1,       0;
          cos(b),	     0,  sin(b)];
      
    % Transform Coordinates of Nodes, Actuators' Tops, Actuators' Bottoms
    try
        for i = 1:node_num
            Nodes.Pos(i,:) = (A*B*Nodes.Pos(i,:)')';
            Actuators.TopPos(i,:) = (A*B*Actuators.TopPos(i,:)')';
            Actuators.BottomPos(i,:) = (A*B*Actuators.BottomPos(i,:)')';
        end
    catch
        disp("Warning: Coordinate System Transformation failed.");
        return;
    end
    a = 0;  b = pi/2;
end
%% About Components
function importData(filepath1, filepath2, filepath3)
    global Nodes;
    global Actuators;
    global Reflectors;
    
    %import data about Crossed Nodes
    opts = detectImportOptions(filepath1);
    opts.VariableTypes = {'string', 'double', 'double', 'double'};
    opts.DataLines = [2, inf];
    opts.SelectedVariableNames = 1;
    Nodes.ID = readmatrix(filepath1, opts);
    opts.SelectedVariableNames = 2:4; 
    Nodes.Pos = readmatrix(filepath1, opts);

    %import data about Actuators
    opts = detectImportOptions(filepath2);
    opts.VariableTypes = {'string', 'double', 'double', 'double', ...
        'double', 'double', 'double'};
    opts.DataLines = [2, inf];
    opts.SelectedVariableNames = 1;
    Actuators.ID = readmatrix(filepath2, opts);
    opts.SelectedVariableNames = 2:4;
    Actuators.BottomPos = readmatrix(filepath2, opts);
    opts.SelectedVariableNames = 5:7;
    Actuators.TopPos = readmatrix(filepath2, opts);
    
    %import data about Reflector elements
    opts = detectImportOptions(filepath3);
    opts.VariableTypes = {'string', 'string', 'string'};
    opts.DataLines = [2, inf];
    opts.SelectedVariableNames = 1:3;
    Reflectors = readmatrix(filepath3, opts);
end
function drawNodes()
    global Nodes;
    plot3(Nodes.Pos(1:end, 1), Nodes.Pos(1:end, 2), Nodes.Pos(1:end, 3),...
        '.k', 'markersize', 8);
end
function drawActuators(node_num)
    global Actuators;
    for i = 1:node_num
        x = [Actuators.TopPos(i, 1), Actuators.BottomPos(i, 1)];
        y = [Actuators.TopPos(i, 2), Actuators.BottomPos(i, 2)];
        z = [Actuators.TopPos(i, 3), Actuators.BottomPos(i, 3)];
        line(x, y, z, 'color', '#333333', 'LineWidth', 5);
    end
end
function drawTiedownCables(node_num)
    global Nodes;
    global Actuators;
    for i = 1:node_num
        x = [Actuators.TopPos(i, 1), Nodes.Pos(i, 1)];
        y = [Actuators.TopPos(i, 2), Nodes.Pos(i, 2)];
        z = [Actuators.TopPos(i, 3), Nodes.Pos(i, 3)];
        line(x, y, z, 'color', '#A2A2A2');
    end
end
function drawReflectors(node_num)
    global Nodes;
    global Reflectors;
    for i = 1:length(Reflectors)
        verts = zeros(3,3);
        for j = 1:node_num
            for k = 1:3
                if strcmp(Reflectors(i, k), Nodes.ID(j))
                    verts(k,:) = Nodes.Pos(j, :);
                end 
            end
        end
        drawTriangle(verts);
    end
    
    function drawTriangle(vertexes)
        for edge = 1:3
            x = [vertexes(edge,1), vertexes(rem(edge, 3) + 1,1)];
            y = [vertexes(edge,2), vertexes(rem(edge, 3) + 1,2)];
            z = [vertexes(edge,3), vertexes(rem(edge, 3) + 1,3)];
            line(x, y, z, 'color', 'k');
        end
    end
end
%% About FAST
%   GETs
function [x, y, z] = getFASTSphCenter()
    x = 0;  y = 0;   z = 0;
end
function getFASTCaliberCenter()
    %TODO
end
function getFASTCaliberCircle()
    %TODO
end
function [X, Y, Z] = getVirtualSphere(R, R_FAST)
    [r, t] = meshgrid(0:R_FAST, 0:0.02:2*pi);
    X = r.*cos(t);
    Y = r.*sin(t);
    Z = - sqrt(R^2 - X.^2 - Y.^2);
end

%   DRAWs
function drawFASTSphCenter()
    plot3(0, 0, 0, 'or', 'markersize', 2);
    text(0, 10, 0, 'C', 'color', 'r');
end
function drawFASTCaliberCenter()
    %TODO
end
function drawFASTCaliberCircle()
    %TODO
end
function drawVirtualSphere(X, Y, Z)
    surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

%% About Light Source
% GETs
function [x, y, z] = getSourcePoint(a, b, r)
    [x, y, z] = sph2cart(a, b, r);
end
function [sol_x, sol_y, sol_z] = getProjectionPoint(sourcePoint, R)
    syms x y z k;
    assume(z<0);
    eqn1 = x == k*sourcePoint(1);
    eqn2 = y == k*sourcePoint(2);
    eqn3 = z == k*sourcePoint(3);
    eqn4 = x^2 + y^2 + z^2 == R^2;
    
    sol = solve([eqn1 eqn2 eqn3 eqn4], [x y z k], 'Real', true);
    sol_x=sol.x;	sol_y=sol.y;	sol_z=sol.z;   
end
% DRAWs
function drawSourcePoint(x, y, z)
    plot3(x, y, z, '.y', 'markersize', 15);
    text(x, y+10, z, 'S', 'color', 'y');
end
function drawProjectionPoint()
    %TODO
end
function drawLightPath(sourcePoint, projPoint)
    % Old Method:
    % line([x -k*x], [y -k*y], [z, -k*z],'color', 'y', 'LineStyle', '--');
    % -----------
    x = [sourcePoint(1), projPoint(1)];
    y = [sourcePoint(2), projPoint(2)];
    z = [sourcePoint(3), projPoint(3)];
    line(x, y, z, 'color', 'y', 'LineStyle', '--');
    
end
%% About Feedback Cabin
% GETs
function [sol_x, sol_y, sol_z] = getCabinCenter(sourcePoint, r)
    syms x y z k
    assume(x<=0 & y<=0 & z<=0);
    eqn1 = x == k*sourcePoint(1);
    eqn2 = y == k*sourcePoint(2);
    eqn3 = z == k*sourcePoint(3);
    eqn4 = x^2 + y^2 + z^2 == r^2;
    sol = solve([eqn1 eqn2 eqn3 eqn4], [x y z k]);
    sol_x=sol.x;	sol_y=sol.y;	sol_z=sol.z;
end
function [X, Y, Z] = getCabinCircle(cabinCenter, sourcePoint, r_cabin)
    [X, Y, Z] = getCircle(cabinCenter, sourcePoint, r_cabin);
end
function [X, Y, Z] = getCabinDisk(cabinCenter, a, b, r_cabin)
    [r, t] = meshgrid(0:r_cabin, 0:0.02:2*pi);
    X = r.*cos(t);
    Y = r.*sin(t);
    Z = zeros(size(X,1), size(X,2)) + cabinCenter(3);
end
% DRAWs
function drawCabinCenter()
    %TODO
end
function drawCabinCircle(X, Y, Z)
    plot3(X, Y, Z,'color', 'b', 'LineWidth', 3);
    text(X(1), Y(1), Z(1)+10, 'Feedback Cabin', 'color', 'b', ...
        'FontSize', 6, 'FontWeight', 'bold');
end
function drawCabinDisk(X, Y, Z)
    % Description:
    %   It appears that the disk is too small to find in the graphic, so we
    % shifted to drawCabinCircle function instead of this one.
    plot3(X, Y, Z, 'color', 'b');
    text(X(1), Y(1), Z(1)+10, 'Feedback Cabin', 'color', 'b', ...
      'FontSize', 6, 'FontWeight', 'bold');
end
%% About Paraboloid
% GETs
function [x, y, z] = getParaCaliberCenter(a, b, R)
    [x, y, z] = sph2cart(a, b, -sqrt(3)/2*R);
end
function [X, Y, Z] = getParaCaliberCircle(paraClbCenter, r_pClb)
    [X, Y, Z] = getCircle(paraClbCenter, paraClbCenter, r_pClb);
end
function [X, Y, Z] = getOptPara(R, F, paraClbCenter)
    % Create Normal Paraboloid
    [r, t] = meshgrid(0:0.5*R, 0:0.02:2*pi);
    X = r.*cos(t);
    Y = r.*sin(t);
    Z = (X.^2 + Y.^2) / (4*F);

    % Translate
    [X, Y, Z] = curveTrans(X, Y, Z, paraClbCenter, R);
    
    % Deserted Steps:
    % 1. Coordinate System Rotation;
    % 2. Sphere Rotation;
end
% DRAWs
function drawParaCaliberCenter()
    %TODO
end
function drawParaCaliberCircle(X, Y, Z)
    plot3(X, Y, Z,'color', '#FFAAFF', 'LineWidth', 3);
    text(X(1), Y(1), Z(1)+10, 'Paraboloid Caliber', 'color', '#FF55FF');
end
function drawOptPara(X, Y, Z)
    surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

%% TOOL FUNCTIONS
function [X, Y, Z] = getCircle(c, n, r)
    theta=(0:2*pi/100:2*pi)'; %theta角从0到2*pi
    a=cross(n,[1 0 0]); %n与i叉乘，求取a向量
    if ~any(a) %如果a为零向量，将n与j叉乘
        a=cross(n,[0 1 0]);
    end
    b=cross(n,a); %求取b向量
    a=a/norm(a); %单位化a向量
    b=b/norm(b); %单位化b向量

    C1=c(1)*ones(size(theta,1),1);
    C2=c(2)*ones(size(theta,1),1);
    C3=c(3)*ones(size(theta,1),1);

    X=C1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%圆上各点的x坐标
    Y=C2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%圆上各点的y坐标
    Z=C3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%圆上各点的z坐标

end
function [X, Y, Z] = curveTrans(X, Y, Z, direction, distance)
%translate the curve
    %Normalize direction vector
    direction = direction / norm(direction);
    
    %Difference Value
    delta = direction * distance;
    
    %Translate Pos
    X = X + delta(1);
    Y = Y + delta(2);
    Z = Z + delta(3);
end