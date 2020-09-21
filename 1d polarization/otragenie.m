clear;
L=520; %������ ������ ������� �� ������� ���������� ���� �����
R0b=10; % ������ �������� ����� ������
N=2^10;
R0=R0b*N/L;
par1=0.01;
lambda=0.63;
k0=2*pi/lambda;
alfa=0.5;
tolshina1=100;
tolshina2=5;
n1=1;
n2=3.5^0.5; %2.131^0.5;
n3=5.8^0.5;
dist=100; %3*lambda; %7.18*lambda;
eps1=1;
eps2=1.5; %2.131;
eps3=5.8;
m=-N/2:1:N/2-1; %����� � �������� ������������ 
h=dist;

i=0;
for ss=-pi*0.5:0.001:pi*0.5
%[m n]=size(xxx);
%for i=1:n
%    ss=pi/2-ugolok(i);
    i=i+1;
    x(i)=ss;
    sinQ(i)=sin(ss);
    sinQ1=sin(ss);
    cosQ1=cos(ss);
    sinQ2=sinQ1*n1/n2;
    cosQ2=(1-sinQ2^2)^0.5;
    sinQ3=sinQ2*n2/n3;
    cosQ3=(1-sinQ3^2)^0.5;
    betta=2*pi/lambda*n2*dist*cosQ2;
    betta1=0;
    betta2=pi/2;
% ����������� ��������� ��� s-����������� (����������������)
r12s=(n1.*cosQ1-n2.*cosQ2)./(n1.*cosQ1+n2.*cosQ2); % ����������� ��������� �� ������� 12
r23s=(n2.*cosQ2-n3.*cosQ3)./(n2.*cosQ2+n3.*cosQ3); % ����������� ��������� �� ������� 23
Rs(i)=(r12s+r23s.*exp(2.*1i.*betta))./(1+r12s.*r23s.*exp(2.*1i.*betta)); % ����������� ����������� ���������
Rs1(i)=(r12s+r23s.*exp(2.*1i.*betta1))./(1+r12s.*r23s.*exp(2.*1i.*betta1)); % ����������� ����������� ���������
Rs2(i)=(r12s+r23s.*exp(2.*1i.*betta2))./(1+r12s.*r23s.*exp(2.*1i.*betta2)); % ����������� ����������� ���������

%����������� ��������� ��� p-����������� (������������)
r12p=(n2.*cosQ1-n1.*cosQ2)./(n2.*cosQ1+n1.*cosQ2); % ����������� ��������� �� ������� 12
r23p=(n3.*cosQ2-n2.*cosQ3)./(n3.*cosQ2+n2.*cosQ3); % ����������� ��������� �� ������� 23
Rp(i)=(r12p+r23p.*exp(2.*1i.*betta))./(1+r12p.*r23p.*exp(2.*1i.*betta)); % ����������� ����������� ���������
Rp1(i)=(r12p+r23p.*exp(2.*1i.*betta1))./(1+r12p.*r23p.*exp(2.*1i.*betta1)); % ����������� ����������� ���������
Rp2(i)=(r12p+r23p.*exp(2.*1i.*betta2))./(1+r12p.*r23p.*exp(2.*1i.*betta2)); % ����������� ����������� ���������

%{
    R12(i)=(n1*cosQ1-n2*cosQ2)/(n1*cosQ1+n2*cosQ2); % ��������� ������������ ������������ ����
    R23(i)=(n2*cosQ2-n3*cosQ3)/(n2*cosQ2+n3*cosQ3); % ��������� ������������ ������������ ����
    RR12(i)=(n1*cosQ1-n2*cosQ2)/(n1*cosQ1+n2*cosQ2); % ��������� ���������������� ������������ ����
    RR23(i)=(n2*cosQ2-n3*cosQ3)/(n2*cosQ2+n3*cosQ3);    
    R(i)=((R12(i)^2+R23(i)^2+2*R12(i)*R23(i)*cos(2*pi/lambda*n2*dist*cosQ2))/(1+R12(i)^2*R23(i)^2+2*R12(i)*R23(i)*cos(2*pi/lambda*n2*dist*cosQ2)))^0.5;
   R(i)=abs((R12(i)+R23(i)*exp(1i*2*betta))/(1+R12(i)*R23(i)*exp(1i*2*betta)));
    RR(i)=((RR12(i)^2+RR23(i)^2+2*RR12(i)*RR23(i)*cos(2*pi/lambda*n2*dist*cosQ2))/(1+RR12(i)^2*RR23(i)^2+2*RR12(i)*RR23(i)*cos(2*pi/lambda*n2*dist*cosQ2)))^0.5;
   RR(i)=abs(RR12(i)+RR23(i)*exp(1i*2*betta)/(1+RR12(i)*RR23(i)*exp(1i*2*betta)));
%}
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���������� ������������� ���������
gam=lambda/L;

gg=(m.* gam).^2;
costheta1= (1-gg).^0.5;
theta=acos(costheta1);
costheta2=(1-(eps1/eps2)*gg).^0.5;

costheta3=(1-(eps1/eps3)*gg).^0.5;

beta=(2*pi/lambda)*n2*h*costheta2;
exb= exp(2*1i*beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����������� ��������������� ��������� �������

r12out=(n1.*costheta1- n2.*costheta2)./ (n1.*costheta1+ n2.*costheta2);

r23out=(n2.*costheta2- n3.*costheta3)./ (n2.*costheta2+ n3.*costheta3);

rout=(r12out+r23out.*exb)./(1+r12out.*r23out.*exb);
rabs=abs(rout);


figure
plot(x,sinQ); 
title('����� x �� -\pi/2 �� \pi/2');
xlabel('���� \alpha');
ylabel('����� ���� sin(\alpha)');

%figure
%plot(x,R12)
%title('���������� ��������� �� ������� ���� 1-2');
%xlabel('���� \alpha');
%ylabel('R_1_2=(n_2*cos(\alpha)-n_1*cos(\beta))/(n_2*cos(\alpha)+n_1*cos(\beta)))');
%figure
%plot(x,R23)
%title('���������� ��������� �� ������� ���� 2-3');
%xlabel('���� \alpha');
%ylabel('R_2_3=(n_3*cos(\beta)-n_2*cos(\gamma))/(n_3*cos(\beta)+n_2*cos(\gamma)))');
figure
plot(x,abs(Rs),'blue')
hold on
plot(x,abs(Rs1),'red')
hold on
plot(x,abs(Rs2),'red')
title('���������� ��������� �� ���������, ����������� ����������� ��������� �������');
%xlabel({'���� \alpha';'';sprintf('���������: n_1=%f, n_2=%f, n_3=%f',n1,n2,n3);sprintf('����� �����= %f ������, ������� ��������� =%f ������',lambda,dist);''});
xlabel({'���� \alpha';'';['���������: n_1=',num2str(n1),', n_2=',num2str(n2),', n_3=',num2str(n3)];['����� ����� \lambda=',num2str(lambda),' ������, ������� ��������� d=',num2str(dist),' ������'];''});
ylabel('R_�_�_�');
%figure
%plot(theta,rabs)
%title('���������� ��������� �� ���������, ����������� ����������� ��������� ������� - ������� ����������');
%xlabel({'���� \alpha';'';sprintf('���������: n_1=%f, n_2=%f, n_3=%f',n1,n2,n3);sprintf('����� �����= %f ������, ������� ��������� =%f ������',lambda,dist);''});
%xlabel({'���� \alpha';'';['���������: n_1=',num2str(n1),', n_2=',num2str(n2),', n_3=',num2str(n3)];['����� ����� \lambda=',num2str(lambda),' ������, ������� ��������� d=',num2str(dist),' ������'];''});
%ylabel('R_�_�_�');
figure
plot(x,abs(Rp))
hold on
plot(x,abs(Rp1),'red')
hold on
plot(x,abs(Rp2),'red')
title('���������� ��������� �� ���������, ����������� ��������������� ��������� �������');
xlabel({'���� \alpha';'';['���������: n_1=',num2str(n1),', n_2=',num2str(n2),', n_3=',num2str(n3)];['����� ����� \lambda=',num2str(lambda),' ������, ������� ��������� d=',num2str(dist),' ������'];''});
ylabel('R_�_�_�');
