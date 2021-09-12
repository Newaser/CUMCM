clc; clear;
%% DATA IMPORT PART
% Import the components
global Nodes;
global Actuators;
global Reflectors;
importComponents('..\Problems\A\附件1.csv', '..\Problems\A\附件2.csv', ...
    '..\Problems\A\附件3.csv');

%% PRETREATMENT PART
% Pre-definition
R = 300;    % Raduis of FAST sphere
R_FAST = 500*0.5;   % Raduis of FAST's caliber
F = 0.466*R;    % Half Focal Length of the Para.
r = R-F;    % Raduis of focal sphere
% alpha_degree = 0;
% beta_degree = 90;
alpha_degree = 36.795;
beta_degree = 78.169;
r_cabin = 1*0.5;
sourceDist = 200;   % Distance from light source to sphere center
RMLim = 0.6;    % Radial Movement Limit
SCERLim = 0.07*0.01;    % Steel Cable Elasticity Rate Limit

% Pre-caculate
node_num = length(Nodes.ID);
alpha = alpha_degree*pi/180;
beta = beta_degree*pi/180;

%% Coordinate System Transformation Part
[alpha, beta] = CoordSysTrans(alpha, beta, node_num);
    

%% Caculation Part
% About Components
    mdfNodes = getModifiableNodes(node_num, R);
    Edges = getEdges(node_num);

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
%% Algorithm Part
% Promblem 2: Radial Approaching Algorithm
% DeltaRhos = RadialApproachingAlgo(mdfNodes, node_num, F, R, Edges, RMLim, SCERLim);
%% Export Part
    exportModifiableNodes(mdfNodes, '.\Exports\modifiable_nodes.xlsx');
    % writematrix(DeltaRhos, '.\Exports\delta_rhos.xlsx');
%% Graphic Plot Part
hold on
% About Components
    drawNodes();
    drawActuators(node_num);
    drawTiedownCables(node_num);
    % drawReflectors(node_num); % Warning: ...
    drawReflectors2(Edges);

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
%% About Algorithm
% About Problem 2:
function delta_rhos = RadialApproachingAlgo(mdfNodes, node_num, F, R, Edges, ...
    RMLim, SCERLim)
% Algorithm ID: #1
    global Nodes;
%     Nodes.ID = mdfNodes.ID;
%     Nodes.Pos = mdfNodes.Pos;
%     node_num = mdfNodes.num;
%     Edges = getEdges(node_num);

    % Define n, m
    n = node_num;
    m = mdfNodes.num;
    
    % Get α β ρ1
    [A, B, rho1] = cart2sph(Nodes.Pos(:,1),Nodes.Pos(:,2),Nodes.Pos(:,3));
    A = A';
    B = B';
    rho1 = rho1';
    
    % Get IM, IE
    IM = mdfNodes.Index;
    IE = Edges.Index;
    
    % Get rhoPara
    rhoPara = zeros(1,m);
    cache_path = '.\Caches\Algorithms\#1_rhoPara.csv';
    cached = true;
    try readmatrix(cache_path);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:textio:textio:FileNotFound')
            cached = false;
        end
    end 
    if cached
        rhoPara = readmatrix(cache_path);
    else
        for i = 1:m
            syms x
            eqn1 =	4*F*R + 4*F*x*sin(B(IM(i))) == ...
                    x^2 * cos(B(IM(i)))^2;
            eqn2 =  x >= 0;
            rhoPara(i) = solve([eqn1 eqn2], x, 'Real', true);
        end
        writematrix(rhoPara, cache_path);         
    end
        
    % Optimal Problem %
    % Optimal Vars:
    d_rho = optimvar('d_rho', n, 'LowerBound',-RMLim,'UpperBound',RMLim);
    
    % ρ2
    rho2 = rho1' + d_rho;
    
    % L1
    L1 = zeros(n,1);
    for idx = 1:n
        a = IE(idx, 1);
        b = IE(idx, 2);
        L1(idx,1) = ...
            sphDist([A(a), B(a), rho1(a)], [A(b), B(b), rho1(b)]);
    end
    
    % L2
    L2 = optimexpr(n);
    for idx = 1:n
        a = IE(idx, 1);
        b = IE(idx, 2);
        L2(idx,1) = ...
            sphDist([A(a), B(a), rho2(a)], [A(b), B(b), rho2(b)]);
    end
        
    % Objective Function
    obj_fun = 0;
    for i = 1:m
        obj_fun = obj_fun + (rho2(IM(i),1) - rhoPara(i))^2;
    end
    
    % Constraints
    SCERUpperLim =  optimconstr(n,1);
    SCERLowerLim =  optimconstr(n,1);
    
    for i = 1:n
        SCERUpperLim(i) =  (L1(i)-L2(i)) / L1(i) <=  SCERLim;
        SCERLowerLim(i) =  (L1(i)-L2(i)) / L1(i) >= -SCERLim;
    end
    
    % Construct the Problem
    prob = optimproblem;
    prob.Objective = obj_fun;
    prob.Constraints.SCERUpperLim = SCERUpperLim;
    prob.Constraints.SCERLowerLim = SCERLowerLim;
    
    % Set Start Point
    x0.d_rho = zeros(n,1);
    
    % Solve
    sol = solve(prob, x0);
    
    %Return value
    delta_rhos = sol.d_rho;
    
%   show(sol);
    disp("OK");
    

%     function lenMatrix = L(rho_i)
%     %Caculate length of edge
%         % Edge Length Column Matrix
%         lenMatrix = zeros(n,1);
%         for idx = 1:n
%             a = IE(idx, 1);
%             b = IE(idx, 2);
%             lenMatrix(idx,1) = ...
%                 sphDist([A(a), B(a), rho_i(a)], [A(b), B(b), rho_i(b)]);
%         end
%     end
% 
%     function sum = obj_fun(rho)
%         sum = 0;
%         for c = 1:m
%             sum = sum + (rho(IM(c),1) - rhoPara(c))^2;
%         end
%     end
end
%% About I/O
%	Input
function importComponents(filepath1, filepath2, filepath3)
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
%   Output
function exportModifiableNodes(mdfNodes, filepath)
    writematrix('节点编号', filepath,'Range','A1', ...
        'WriteMode', 'overwritesheet');
    writematrix('X坐标(米)', filepath,'Range','B1');
    writematrix('Y坐标(米)', filepath,'Range','C1');
    writematrix('Z坐标(米)', filepath,'Range','D1');
    writematrix('θ坐标(rad)', filepath,'Range','E1');
    writematrix('γ坐标(rad)', filepath,'Range','F1');
    writematrix('r坐标(米)', filepath,'Range','G1');
    writematrix(mdfNodes.ID',filepath,'Range','A2');
    writematrix(mdfNodes.Pos,filepath,'Range','B2');
    [theta, gamma, r] = cart2sph(mdfNodes.Pos(:,1),mdfNodes.Pos(:,2),mdfNodes.Pos(:,3));
    writematrix([theta, gamma, r],filepath,'Range','E2');
    
end
%% About Coordinate System
function [a, b] = CoordSysTrans(a, b, node_num)
    global Nodes;
    global Actuators;
    log_path = '.\Caches\last_ab.txt';
    
    % Reset some data
    resetData();
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
    writematrix([a b], log_path);
    a = 0;  b = pi/2;
    
    function resetData()
        % Delete caches of Algorithms
        last_ab = readmatrix(log_path);
        if (a-last_ab(1))^2 + (b-last_ab(2))^2 > 10^(-5)
            rmdir('.\Caches\Algorithms', 's');
            mkdir('.\Caches\Algorithms');
        end       
    end
end
%% About Components
function movedNodes = moveNodes(node_num, IDs, destPos)
    global Nodes;
    
    move_num = length(IDs);
    movedNodes.ID = Nodes.ID;
    movedNodes.Pos = Nodes.Pos;
    for i = 1:move_num
        for j = node_num
            if strcmp(IDs(i), movedNodes.ID(i))
                movedNodes.Pos(i,:) = destPos(i,:);
            end
        end
    end
end
function mdfNodes = getModifiableNodes(node_num, R)
    global Nodes;
    
    mdfNodes.ID = strings(node_num);
    mdfNodes.Pos = zeros(node_num,3);
    mdfNodes.Index = zeros(node_num);
    mdfNodes.num = 0;
    
    for i = 1:node_num
        if Nodes.Pos(i,3) <= -sqrt(3)/2 * R
            mdfNodes.num = mdfNodes.num + 1;
            mdfNodes.ID(mdfNodes.num) = Nodes.ID(i);
            mdfNodes.Pos(mdfNodes.num,:) = Nodes.Pos(i,:);
            mdfNodes.Index(mdfNodes.num) = i;
        end
    end
    
    mdfNodes.ID = mdfNodes.ID(1:mdfNodes.num);
    mdfNodes.Pos = mdfNodes.Pos(1:mdfNodes.num,:);
    mdfNodes.Index = mdfNodes.Index(1:mdfNodes.num);
end
function Edges = getEdges(node_num)
    global Nodes;
    global Reflectors;
    
    % Init
    Edges.IDs = strings(node_num*3,2);
    Edges.Index = zeros(node_num*3,2);
    Edges.num = 0;
    cached = true;
    
    try readmatrix('.\Caches\edges.xlsx');
    catch ME
        if strcmp(ME.identifier, 'MATLAB:textio:textio:FileNotFound')
            cached = false;
        end
    end
    
    if cached
        Edges.num = readmatrix('.\Caches\edges.xlsx', 'Range', 'A1:A1');
        Edges.Index = readmatrix('.\Caches\edges.xlsx', 'Range', 'B1');
        
        Edges.IDs = strings(Edges.num,2);
        for i = 1:Edges.num
            for j = [1 2]
                Edges.IDs(i, j) = Nodes.ID(Edges.Index(i,j));
            end
        end
    else
        for i = 1:size(Reflectors,1)
            for j = 1:3
                vertex1 = Reflectors(i,j);
                vertex2 = Reflectors(i, rem(j, 3) + 1);
                if ~existsEdge([vertex1, vertex2])
                    Edges.num = Edges.num + 1;
                    Edges.IDs(Edges.num,:) = [vertex1, vertex2];
                    Edges.Index(Edges.num,:) = ...
                        [nodeIndex(vertex1), nodeIndex(vertex2)];
                end
            end
        end

        Edges.IDs = Edges.IDs(1:Edges.num,:);
        Edges.Index = Edges.Index(1:Edges.num,:);
        writematrix(Edges.num,'.\Caches\edges.xlsx', 'Range', 'A1:A1');
        writematrix(Edges.Index,'.\Caches\edges.xlsx', 'Range', 'B1');
        
    end
    
    function idx = nodeIndex(nodeID)
        for c = 1:node_num
            if strcmp(nodeID, Nodes.ID(c))
                idx = c;
                return;
            end
        end
        idx = 0;
    end
    function yes = existsEdge(e)
        for a = 1:Edges.num
            if edgeEqual(e, edge(a))
                yes = true;
                return;
            end
        end
        yes = false;
    end
    function yes = edgeEqual(e1, e2)
        if ( strcmp(e1(1),e2(1)) && strcmp(e1(2),e2(2)) ) || ...
           ( strcmp(e1(1),e2(2)) && strcmp(e1(2),e2(1)) )
            yes = true;
        else
            yes = false;
        end
    end
    function e = edge(i)
        e = Edges.IDs(i,:);
    end
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
function drawReflectors2(Edges)
    global Nodes;
    for i = 1:Edges.num
        x = [Nodes.Pos(Edges.Index(i,1),1), Nodes.Pos(Edges.Index(i,2),1)];
        y = [Nodes.Pos(Edges.Index(i,1),2), Nodes.Pos(Edges.Index(i,2),2)];
        z = [Nodes.Pos(Edges.Index(i,1),3), Nodes.Pos(Edges.Index(i,2),3)];
        line(x, y, z, 'color', 'k', 'LineWidth', 0.5);
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
    % Z = (X.^2 + Y.^2 - 1.6794*10^5) / (2*279.9038) + R;
    
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
function dist = sphDist(p1, p2)
    [x1, y1, z1] = sph2cart(p1(1), p1(2), p1(3));
    [x2, y2, z2] = sph2cart(p2(1), p2(2), p2(3));
    dist = norm([x1, y1, z1]-[x2, y2, z2]);
end