density = 1065;
Told = 40;
[nuOld,DTnu] = kinematic_viscosity([],density,Told,1);

nuNew = nan(2,1);
del = 1e-7;
nuNew(1) = kinematic_viscosity([],density,Told-del,1);
nuNew(2) = kinematic_viscosity([],density,Told+del,1);
FD_nu = 0.5*(nuNew(2) - nuNew(1))/del;

fprintf('absolute difference = %g \n', max(abs(DTnu - FD_nu)))