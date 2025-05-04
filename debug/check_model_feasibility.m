

f =  sdpvar(length(paths),1);

% Constraints
cnstr = [X_paths_G_reduced*f >= X_copy_od(arc_mask)-epsl_od/2,...
        X_paths_G_reduced*f <= X_copy_od(arc_mask)+epsl_od/2, 
        f>=0];

model.A = sparse([X_paths_G_reduced;X_paths_G_reduced;1]);
model.obj = [0];
% model.rhs = [X_paths_G_reduced.*X_copy_od(arc_mask);%-X_paths_G_reduced*epsl_od;
%              X_paths_G_reduced.*X_copy_od(arc_mask);%+X_paths_G_reduced*epsl_od;
%              0];
model.rhs = [X_copy_od(arc_mask)-epsl_od*10;
             X_copy_od(arc_mask)+epsl_od*10;
             0];
model.sense = [repmat('>', size(X_paths_G_reduced));repmat('<', size(X_paths_G_reduced)); '>'];

