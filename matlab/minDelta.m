function [minDelta,C] = SDPnormtest
tic
N    = 8;
din  = 2;
dout = 2;

rho_i = zeros(dout,dout,N);
sigma_i = zeros(din,din,N);

for i=1:N
    rho_i(:,:,i)=RandomDensityMatrix(dout)
    sigma_i(:,:,i)=RandomDensityMatrix(din)
end

yalmip('clear')

C = sdpvar(din*dout,din*dout,'hermitian','complex');

F = [C>=0,PartialTrace(C,2)==eye(din)];

minDelta = 0;
for i=1:N
    minDelta = minDelta + norm(rho_i(:,:,i)-PartialTrace(kron(transpose(sigma_i(:,:,i)),eye(dout))*C,[1],[din dout]),1);
end
minDelta=minDelta/N;

ops = sdpsettings( ...
        'solver'       , 'sedumi', ...
        'sedumi.eps'   , 1e-8,       ... % tolerances
        'sedumi.maxiter', 300 , ...
        'sedumi.stepdif', 0   , ... % show every interior point step
        'verbose'      , 2,   ...
        'debug'        , 1,   ...   % keeps temporary files
        'savesolverinput',1, ...
        'savesolveroutput',1);


% solution = solvesdp(F,minDelta,sdpsettings('verbose',1,'cachesolvers',1));
%solution = solvesdp(F,minDelta,sdpsettings('verbose',1,'cachesolvers',1,'solver','sedumi'));
solution = solvesdp(F, minDelta, sdpsettings(ops))

C = double(C);
minDelta = double(minDelta);
total_time = toc
end
