% try readmatrix('not_exists.txt');
% catch ME
%     if strcmp(ME.identifier, 'MATLAB:textio:textio:FileNotFound')
%     end
% end
% %writematrix([1 2],'here.txt');

A = [1 2 4;6 7 8];
disp(A);
disp(A');

B = [1;2;4];
disp(B);
disp(B');