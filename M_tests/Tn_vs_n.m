ftsz = 30;
normp = [4,8,12,16,20,24,32,64,128,256,512];
area = 0.15*0.2;
triNpt1D = 4; % for line integration in regular FEM,poly IGFEM and NURBS IGFEM
triNpt2D = 61; % for element integration in regular FEM and poly IGFEM
gauss2.line = gauss_points_and_weights(true,triNpt1D,1,'combined');
gauss2.elem = gauss_points_and_weights(true,triNpt2D ,2,'combined');
Tmax = max(UUR2);
nNormp = numel(normp);
pnormT = nan(nNormp,1);
pnormT1 = nan(nNormp,1);
pnormTA = nan(nNormp,1);
for i = 1:nNormp
    pnormT(i) = field_p_norm(mesh.elem,mesh.node.coords,UUR,...
                             normp(i),Tmax,gauss2);
    pnormT1(i) = field_p_norm(mesh.elem,mesh.node.coords,UUR,...
                             normp(i),1,gauss2); 
    pnormTA(i) = pnormT(i)/area^(1/normp(i)); 
end



figure
plot(normp,pnormT,'color','r','linestyle','-','marker','o','linewidth',2)
hold on
plot(normp,pnormT1,'color','b','linestyle','-.','marker','>','linewidth',2)
plot(normp,pnormTA,'color','g','linestyle','-','marker','^','linewidth',2)
plot([min(normp),max(normp)],[Tmax,Tmax],'color','k','linestyle','--','linewidth',2)
xlabel('n','fontsize',ftsz)
ylabel('T','fontsize',ftsz)
legend('T_n','T_{max}')
set(gca,'fontsize',ftsz)