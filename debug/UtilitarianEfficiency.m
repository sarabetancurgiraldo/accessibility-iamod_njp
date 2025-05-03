function UtilitarianEfficiency(nCar,G,B,nArcs,D,str_save)
%% Define parameters
arcsCar                 = find(G.Edges.Type == 1);
nCarArcs                = sum(G.Edges.Type == 1);
Bcar                    = B(:,arcsCar);
t                       = G.Edges.Weight; % in hours

%% Define variables
X                       = sdpvar(nArcs, size(D,2), 'full');
xR                      = sdpvar(nCarArcs, 1, 'full');

%% Optimization

% Define constraints
Cons                    = [B*X                                   == D;
                           Bcar*(sum(X(arcsCar,:),2)+xR)         == 0;
                           t(arcsCar)'*(sum(X(arcsCar,:),2)+xR)  <= nCar;
                           X                                     >= 0;
                           xR                                    >= 0];

% Define objective
Obj                     = t'*X*ones(size(D,2), 1); 

options                 = sdpsettings('verbose', 1, ...
                                      'solver', 'gurobi', ...
                                      'showprogress', 1);

%% Save on object
sol_utilEff             = optimize(Cons, Obj, options);
sol_utilEff.X           = value(X);
xR                      = value(xR);
xRfull                  = zeros(nArcs,1);
xRfull(arcsCar)         = xR;
sol_utilEff.xR          = xRfull;

% save data
save(str_save,"sol_utilEff");
yalmip('clear');

end