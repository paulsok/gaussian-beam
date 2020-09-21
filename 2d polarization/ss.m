hold all
plot(mkm,El(512,:),'LineWidth',3,'Color','green');
plot(mkm,El(512,:),'LineWidth',3,'Color','red');
title('Сечение по центру пучка (до и после линзы)','fontsize',12)
xlabel({'Пространственная координата,мкм'},'fontsize',12)
xlim([llim rlim])
ylabel('Угол эллиптичности  ^0','fontsize',12)

hold all
plot(mkm,Az(512,:),'LineWidth',3,'Color','green');
plot(mkm,Az(512,:),'LineWidth',3,'Color','red');
title('Сечение по центру пучка (до и после линзы)','fontsize',12)
xlabel({'Пространственная координата,мкм'},'fontsize',12)
ylabel('Азимут угла ^0','fontsize',12);
xlim([llim rlim]) 



