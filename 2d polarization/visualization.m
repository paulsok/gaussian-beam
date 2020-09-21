%% ������������.
%% >����� (����������� & �������������)------------------------------------------------------------%
subplot(2,3,1);
plot(mkm,El(512,:),'LineWidth',3,'Color','red');
title('������� �� ������ �����','fontsize',12)
xlabel({'���������������� ����������,���'},'fontsize',12)
xlim([llim rlim])
ylabel('���� �������������  ^0','fontsize',12)

subplot(2,3,2);
plot(mkm,Az(512,:),'LineWidth',3,'Color','green');
title('������� �� ������ �����','fontsize',12)
xlabel({'���������������� ����������,���'},'fontsize',12)
ylabel('������ ���� ^0','fontsize',12);
xlim([llim rlim]) 

subplot(2,3,3);
plot(mkm,Int4(512,:),'LineWidth',3);
title('������� �� ������ �����','fontsize',12)
xlabel({'���������������� ����������,���'},'fontsize',12)
xlim([llim rlim]) 
ylabel({'�������������,��/�^2'},'fontsize',12);

subplot(2,3,4);
imagesc(El)
title('���������� ������� �����,���� �������������^0','fontsize',12)
axis off
colorbar('location','SouthOutside')

subplot(2,3,5);
imagesc(Az)
title('���������� ������� �����, ������ ����^0','fontsize',12)
axis off
colorbar('location','SouthOutside')

subplot(2,3,6);
surf(X(1:3:N,1:3:N),Y(1:3:N,1:3:N),Int4(1:3:N,1:3:N)/max(max(Int4)))
colormap(gray)
view(0, 90);
title('������������� I / I_m_a_x','fontsize',12)
shading interp 
colorbar('location','SouthOutside')
lighting phong
xlabel({'X, ���'},'fontsize',12)
ylabel({'Y, ���'},'fontsize',12)
xlim([llim rlim]) 
ylim([llim rlim])


%% >������������ ��������� ------------------------------------------------------------------------%
%{
figure
plot(Q,abs(Rp(512,:)),Q,abs(Rs(512,:)),'LineWidth',2)
xlabel({'���� �������, �������'},'fontsize',12)
ylabel({'������������ ��������� |R_s|, |R_p|'},'fontsize',12);
legend('|R_p|','|R_s|','Location','Best','fontsize',12)

figure
imagesc(abs(Rs));figure(gcf);
title('����������� ��������� |R_s|','fontsize',14)
xlabel('X','fontsize',12)
ylabel('Y','fontsize',12);
colorbar
figure
imagesc(abs(Rp));figure(gcf);
title('����������� ��������� |R_p|','fontsize',14)
xlabel('X','fontsize',12)
ylabel('Y','fontsize',12);
colorbar
%}
%% >������������ ����������� ----------------------------------------------------------------------%
%{
figure
hold all
plot(Q,abs(Ts(512,:))*100,Q,abs(Tp(512,:))*100,'LineWidth',2)
xlabel({'���� �������, �������'},'fontsize',12)
ylabel({'������������ ����������� |T_s|, |T_p|'},'fontsize',12)
legend('|T_s|','|T_p|','Location','Best')

figure
imagesc(abs(Ts));figure(gcf);
title('����������� ����������� |T_s|','fontsize',14)
xlabel('X','fontsize',12)
ylabel('Y','fontsize',12);
colorbar
figure
imagesc(abs(Tp));figure(gcf);
title('����������� ����������� |T_p|','fontsize',14)
xlabel('X','fontsize',12)
ylabel('Y','fontsize',12);
colorbar
%}
