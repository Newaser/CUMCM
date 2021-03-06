clc; clear;
%% DATA IMPORT PART
% Import the components
global Nodes;
global Actuators;
global Reflectors;
importData('..\Problems\A\附件1.csv', '..\Problems\A\附件2.csv', ...
    '..\Problems\A\附件3.csv');

%% PRETREATMENT PART
% pre-definition
R = 300;
F = 0.466*R;
r = R-F;
alpha_degree = 0;
beta_degree = 90;
% alpha_degree = 36.795;
% beta_degree = 78.169;
r_cabin = 0.5;

% pre-caculate
node_num = length(Nodes.ID);
alpha = alpha_degree*pi/180;
beta = beta_degree*pi/180;

%% Caculation Part
% About FAST
    getFASTSphCenter();
    getFASTCaliberCenter();
    getFASTCaliberCircle();

% About Light Source
    getSourcePoint();
    getProjectionPoint();

% About Feedback Cabin
    getCabinCenter();
    getCabinDisk();

% About Paraboloid
    getParaCaliberCenter();
    getParaCaliberCircle();
    getOptPara();


%% Graphic Plot Part
hold on
% About Components
    drawNodes();
    drawActuators();
    drawTiedownCables();
    drawReflectors();

% About FAST
    drawFASTSphCenter();
    % drawFASTCaliberCenter();
    drawFASTCaliberCircle();

% About Light Source
    drawSourcePoint();
    % drawProjectionPoint();
    drawLightPath();

% About Feedback Cabin
    % drawCabinCenter();
    drawCabinDisk();
    
% About Paraboloid
    drawParaCaliberCenter();
    drawParaCaliberCircle();
    drawOptPara();

%draw all the nodes
plot3(Nodes.Pos(1:end, 1), Nodes.Pos(1:end, 2), Nodes.Pos(1:end, 3), ...
    '.k', 'markersize', 8);

%draw the center point
plot3(0, 0, 0, 'or', 'markersize', 2);
text(0, 10, 0, 'C', 'color', 'r');

%draw the light source
[x, y, z] = drawSource(alpha, beta, 200, 1.5);
sourceCenter = [x, y, z];

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
% for i = 1:length(Reflectors)
%     vertexes = zeros(3,3);
%     for j = 1:node_num
%         for k = 1:3
%             if strcmp(Reflectors(i, k), Nodes.ID(j))
%                 vertexes(k,:) = Nodes.Pos(j, :);
%             end 
%         end
%     end
%     drawTriangle(vertexes);
% end

%draw focal hemisphere
% drawFocalHemi(r);

%draw feedback cabin
drawFeedbackCabin(sourceCenter, r, r_cabin);

%draw the caliber circle
 [x, y, z] = drawCaliber(alpha, beta, R);
 caliberCenter = [x, y, z];

%draw the paraboliod
drawPara(R, F, alpha, beta, caliberCenter);
%% 

%% 
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
    [x, y, z] = sph2cart(a, b, r);
    plot3(x, y, z, '.y', 'markersize', 15);
    text(x, y+10, z, 'S', 'color', 'y');
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
    Sol = solve([eqn1 eqn2 eqn3 eqn4], [x y z k]);
    cabinCenter = [Sol.x, Sol.y, Sol.z];
    [X, Y, Z] = getCircle(S, r_cabin, cabinCenter);
    plot3(X, Y, Z,'color', 'b', 'LineWidth', 3);
    text(Sol.x, Sol.y+10, Sol.z, 'Feedback Cabin', 'color', 'b', ...
        'FontSize', 6, 'FontWeight', 'bold');
end

function [x, y, z] = drawCaliber(a, b, R)
    [x, y, z] = sph2cart(a, b, -sqrt(3)/2*R);
    [X, Y, Z] = getCircle([x, y, z], R/2, [x, y, z]);
    plot3(X, Y, Z,'color', '#FFAAFF', 'LineWidth', 3);
    text(X(1), Y(1), Z(1)+10, 'Paraboloid Caliber', 'color', '#FF55FF');
end

function drawPara(R, F, a, b, caliberCenter)
%draw the Paraboloid
    %Create Scatter set
    [r, t] = meshgrid(0:0.5*R, 0:0.02:2*pi);
    X = r.*cos(t);
    Y = r.*sin(t);
    Z = (X.^2 + Y.^2) / (4*F);

    %Translate
    [X, Y, Z] = curveTrans(X, Y, Z, caliberCenter, R);
    
    %Coordinate System Rotation
    
    
    %Draw Normal Para.
    Para = surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
%     %Rotate
%     sphRotate(Para, a, b);
end

function [X, Y, Z] = getCircle(n, r, c)
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

function sphRotate(fig, a, b)
% rotate the fig with azimuth a and elevation b

    %Create 2 sample vectors
    [x1,y1,z1] = sph2cart(a,0,1);
    [x2,y2,z2] = sph2cart(a,b,1);
    vec1 = [x1,y1,z1];
    vec2 = [x2,y2,z2];
    
%     %Obtain the Rotation Angle
%     theta = acos(vec1*vec2' / (norm(vec1)*norm(vec2)) );
%     theta_degree = theta*180/pi;
    
    %Obtain the Rotation Axis
    syms x y;
    eqn1 = [x,y,1]*vec1' == 0;
    eqn2 = [x,y,1]*vec2' == 0; 
    sol = solve([eqn1, eqn2], [x, y]);
    direction = [sol.x,sol.y,1];
    
    %Rotation
    rotate(fig, direction, b);
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
%% 

