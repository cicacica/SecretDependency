function [minDelta,C] = SDPnormtest
tic
N    = 8;
din  = 2;
dout = 2;

data = load("experimenttomo.mat")
% this loaded as [8 2 2] matrices
rtomo = data.Expression1.rtomo;
rclean = data.Expression1.rclean;

rho_i = zeros(dout,dout,N);
sigma_i = zeros(din,din,N);

for i = 1:N
    rho_i(:,:,i) = rtomo(i, :, :);
    sigma_i(:,:,i) = rclean(i, :, :);
end


%test if the matrices are assigned correctly
%for i = 1:N
%    fprintf('\n--- Pair %d ---\n', i);
%    fprintf('rho_%d:\n', i);
%    disp(rho_i(:,:,i));
%    fprintf('sigma_%d:\n', i);
%    disp(sigma_i(:,:,i));
%end


yalmip('clear')

C = sdpvar(din*dout,din*dout,'hermitian','complex');

F = [C>=0,PartialTrace(C,2)==eye(din)];

minDelta = 0;
for i=1:N
    minDelta = minDelta + norm(rho_i(:,:,i)-PartialTrace(kron(transpose(sigma_i(:,:,i)),eye(dout))*C,[1],[din dout]),1);
end
minDelta=minDelta/N;

%% solution = solvesdp(F,minDelta,sdpsettings('verbose',1,'cachesolvers',1));
solution = solvesdp(F,minDelta,sdpsettings('verbose',1,'cachesolvers',1,'solver','sedumi'));

C = double(C);
minDelta = double(minDelta);
total_time = toc
end
