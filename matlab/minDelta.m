function [minDelta,C] = SDPnormtest
tic
N    = 8;
din  = 2;
dout = 3;

rho_i = zeros(dout,dout,N);
sigma_i = zeros(din,din,N);

for i=1:N
    rho_i(:,:,i)=RandomDensityMatrix(dout);
    sigma_i(:,:,i)=RandomDensityMatrix(din);
end

yalmip('clear')

C = sdpvar(din*dout,din*dout,'hermitian','complex');

F = [C>=0,PartialTrace(C,2)==eye(din)];

minDelta = 0;
for i=1:N
    minDelta = minDelta + norm(rho_i(:,:,i)-PartialTrace(kron(transpose(sigma_i(:,:,i)),eye(dout))*C,[1],[din dout]),1);
end
minDelta=minDelta/N;

% solution = solvesdp(F,minDelta,sdpsettings('verbose',1,'cachesolvers',1));
solution = solvesdp(F,minDelta,sdpsettings('verbose',1,'cachesolvers',1,'solver','sedumi'));

C = double(C);
minDelta = double(minDelta);
total_time = toc
end