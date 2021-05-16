pts = [0,1;
       0,1];
pts = rand(2,2);
%pts = [4.217612826262750e-01     7.922073295595544e-01;
%       9.157355251890671e-01     9.594924263929030e-01];
figure
plot(pts(1,:),pts(2,:),'r-')
vec = pts(:,2) - pts(:,1);
vecLen = norm(vec);
unitVec = vec/vecLen;
vecLenCube = vecLen^3;

DtDpts = zeros(2,4);
DtDpts(1,1) = -(pts(2,1) - pts(2,2))^2;
DtDpts(2,1) = (pts(2,1) - pts(2,2))*(pts(1,1) - pts(1,2));
DtDpts(1,2) = DtDpts(2,1);
DtDpts(2,2) = -(pts(1,1) - pts(1,2))^2;
DtDpts(1,3) = -DtDpts(1,1);
DtDpts(2,3) = -DtDpts(2,1);
DtDpts(1,4) = -DtDpts(2,1);
DtDpts(2,4) = -DtDpts(2,2);
DtDpts = DtDpts/vecLenCube;
disp('DtDpts = ')
disp(DtDpts)

del = 1e-8;
FD_DtDpts = nan(size(DtDpts));
for i = 1:numel(pts)
    newPts = pts;
    newPts(i) = pts(i) + del;
    newVec = newPts(:,2) - newPts(:,1);
    newUnitVec = newVec/norm(newVec);
    FD_DtDpts(:,i) = (newUnitVec-unitVec)/del;
end

disp('FD_DtDpts = ')
disp(FD_DtDpts)

fprintf('max abs diff = %g \n', max(max(abs(FD_DtDpts - DtDpts))))