%{
ave = average_temp(mesh.elem,mesh.node.coords,UUR2);
fprintf('average temp = %g \n', ave);
norm1 = field_p_norm(mesh.elem,mesh.node.coords,UUR,1,gauss);
fprintf('1-norm = %g \n', norm1);
estnorm1 = field_p_norm(mesh.elem,mesh.node.coords,UUR,1,[]);
fprintf('est 1-norm = %g \n', norm1);
%}
%{
p = 120;
normp = field_p_norm(mesh.elem,mesh.node.coords,UUR,p,gauss);
fprintf('%i-norm = %g \n', p,normp);
estnormp = field_p_norm(mesh.elem,mesh.node.coords,UUR,p,[]);
fprintf('est %i-norm = %g \n',p, normp);
fprintf('relative difference = %g \n',abs(1-estnormp/normp))
%}
p = [1,2,4,8,16,32,64,128,256,512];
normp = zeros(numel(p),1);
estnormp = zeros(numel(p),1);

time1 = 0;
time2 = 0;
for i = 1:numel(p)
    tic;
    normp(i) = field_p_norm(mesh.elem,mesh.node.coords,UUR,p(i),gauss);
    time1 = time1 + toc;
    tic
    estnormp(i) = field_p_norm(mesh.elem,mesh.node.coords,UUR,p(i),[]);
    time2 = time2 + toc;
end
fprintf('time1 = %g, time2 = %g\n',time1,time2)

ftsize = 20;
figure
plot(p,normp,'k-',p,estnormp,'b--','linewidth',2)
hold on
maxField = ones(numel(p),1)*max(UUR2);
plot(p,maxField,'r-.','linewidth',2)
xlabel('p','fontsize',ftsize)
ylabel('p-norm','fontsize',ftsize)
set(gca(),'fontsize',ftsize)
legend('gauss','area','max')


reldiff = abs(1-estnormp./normp);
figure
plot(p,reldiff,'k-','linewidth',2)
xlabel('p','fontsize',ftsize)
ylabel('rel diff','fontsize',ftsize)
set(gca(),'fontsize',ftsize)