hold all
plot(mkm,El(512,:),'LineWidth',3,'Color','green');
plot(mkm,El(512,:),'LineWidth',3,'Color','red');
title('������� �� ������ ����� (�� � ����� �����)','fontsize',12)
xlabel({'���������������� ����������,���'},'fontsize',12)
xlim([llim rlim])
ylabel('���� �������������  ^0','fontsize',12)

hold all
plot(mkm,Az(512,:),'LineWidth',3,'Color','green');
plot(mkm,Az(512,:),'LineWidth',3,'Color','red');
title('������� �� ������ ����� (�� � ����� �����)','fontsize',12)
xlabel({'���������������� ����������,���'},'fontsize',12)
ylabel('������ ���� ^0','fontsize',12);
xlim([llim rlim]) 



