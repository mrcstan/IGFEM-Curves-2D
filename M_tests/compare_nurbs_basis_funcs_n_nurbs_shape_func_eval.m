clear
nurbsInd{1} = [2,2];
nurbsInd{2} = [3,2];
nurbsInd{3} = [3,1];
nurbsInd{4} = [3,3];
nurbsInd{5} = [1,2];
nurbsInd{6} = [2,1];
nurbsInd{7} = [1,1];
nurbsInd{8} = [1,3];
nurbsInd{9} = [2,3];
%{
p(:,:,1) = [9.999999999999995e-03     9.550258919995010e-03     1.018175337102023e-02;
            0     4.343084059617441e-03     8.295515456129925e-03;
            0                         0                         0;
            1.000000000000000e+00     9.550258919995014e-01     9.238634413575078e-01];
        
p(:,:,2) = [1.000000000000000e-02     1.000000000000000e-02     1.000000000000000e-02;
            1.000000000000000e-02     1.000000000000000e-02     1.000000000000000e-02;
                         0                         0                         0;
            1.000000000000000e+00     1.000000000000000e+00     1.000000000000000e+00];

knots = {[0,0,0,1,1,1],[0,0,1,1]};
%}
p(:,:,1) = [0,0;
            0.5,1;
            1,0]';
p(:,:,3) = [0,3;
            0.5,4;
            1,3]';
p(:,:,2) = [0,2;
            0.5,2.2;
            1,2]';

knots = {[0,0,0,1,1,1],[0,0,0,1,1,1]};
nurbsSurf = nrbmak(p,knots);
dnurbsSurf = nrbderiv(nurbsSurf);
[~,J] = nrbdeval(nurbsSurf,dnurbsSurf,[0.2,0.5;0.3,0.1]);

nurbs_plot(nurbsSurf,[50,50]);
ftsize = 20;
xlabel('x','fontsize',ftsize)
ylabel('y','fontsize',ftsize)
zlabel('z','fontsize',ftsize)
set(gca(),'fontsize',ftsize)

%param = {0.5,0.5};
%param2 = rand(2,1)
param2 = [0.1;0];
[B,ind] = nrbbasisfun(param2,nurbsSurf);
[dBu,dBv,ind] = nrbbasisfunder(param2,nurbsSurf);

[nurbsShape,dNurbsShape] = create_nurbs_enrichment_function(nurbsSurf,nurbsInd);
[N,DN] = nurbs_shape_func_eval(nurbsShape,dNurbsShape,param2);

ind2 = zeros(numel(nurbsInd),1);
for i = 1:numel(nurbsInd)
    N(i) = N(i)*nurbsSurf.coefs(4,nurbsInd{i}(1),nurbsInd{i}(2));
    DN(i,:) = DN(i,:)*nurbsSurf.coefs(4,nurbsInd{i}(1),nurbsInd{i}(2));
    ind2(i)= sub2ind(nurbsSurf.number,nurbsInd{i}(1),nurbsInd{i}(2));
end


fprintf('max error of basis functions = %g \n',max(abs(N-B(ind2)')));
fprintf('max relative error of derivative of basis functions in u direction = %g \n',max(abs(DN(:,1)-dBu(ind2)')));
fprintf('max relative error of derivative of basis functions in v direction = %g \n',max(abs(DN(:,2)-dBv(ind2)')));


