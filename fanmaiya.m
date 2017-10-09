clear all;
close all;
clc;


%--------------------------------------------------------------------------
GHz = 1e9; MHz = 1e6;  KHz = 1e3; Hz = 1;
us = 1e-6;  ms = 1e-3;  ns = 1e-9; ps =1e-12;
km = 1e3;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


t_s = -tau_s/2 : 1/fs : tau_s/2 - 1/fs;
k_s = b/tau_s;

weizhi1 = 5000;
weizhi2 = 6000;
piancha = 300;
x = zeros(1,fft_changdu);
% fft_changdu = length(x);
S=exp(1j*(pi*k_s*t_s.^2)); %���Ե�Ƶ�ź�
x_ppxs=fft(S,fft_changdu);%ƥ��ϵ��
flag = 1;

x(1,weizhi1) = 1;
x(1,weizhi2) = 1;% �ӵڶ���Ŀ��
x = x + rand(1,fft_changdu)/10;% ���������
if flag~=1;
    figure,plot(real(x));
    hold on,plot(imag(x),'r');
    xlim([weizhi-piancha,weizhi+piancha]);
end

s1 = fft(x, fft_changdu);% ����ѹ
s2 = conj(s1.*x_ppxs);       % fft֮������ݽ��е�˲���
y = conj(fft((s2),fft_changdu));  % ����ѹ��Ҫ�������ι���Ż��Ϊ��ȷλ�ã���һ����s2���ڶ�����y
% y = fft(s2,length(x));  
if flag~=1;
    figure,plot(real(y));
    hold on,plot(imag(y),'r');
    xlim([weizhi-piancha,weizhi+piancha]);
end

s3 = fft(y, fft_changdu);% ��ѹ
s4 = conj(s3).*x_ppxs;       % fft֮������ݽ��й������
z = fft(s4,fft_changdu);  
if flag~=1;
    figure,plot(real(z));
    hold on,plot(imag(z),'r');
    xlim([weizhi-piancha,weizhi+piancha]);
end

if flag==1
    figure,
    subplot(3,1,1),plot(real(x)),
    hold on,plot(imag(x),'r'),
%     xlim([weizhi-piancha,weizhi+piancha]),
    title('ģ��Ŀ��λ�ã���ɫʵ����ɫ�鲿��'),
    hold off;
    subplot(3,1,2),plot(real(y)),
    hold on,plot(imag(y),'r'),
%     xlim([weizhi-piancha,weizhi+piancha]),
    title('����ѹ��������ɫʵ����ɫ�鲿��'),
    hold off;
    subplot(3,1,3),plot(real(z)),
    hold on,plot(imag(z),'r'),
%     xlim([weizhi-piancha,weizhi+piancha]),
  	title('��ѹ��������ɫʵ����ɫ�鲿��'),
	hold off;
end

return

fp1=fopen('test_zhu.dat','wt');
for n = 1 : changdu
    fprintf(fp1,'%22.14e\n',real(test_zhu(n)));
    fprintf(fp1,'%22.14e\n',imag(test_zhu(n)));
end
fclose(fp1);
fp1=fopen('test_fu1.dat','wt');
for n = 1 : changdu
    fprintf(fp1,'%22.14e\n',real(test_fu1(n)));
    fprintf(fp1,'%22.14e\n',imag(test_fu1(n)));
end
fclose(fp1);
fp1=fopen('test_fu2.dat','wt');
for n = 1 : changdu
    fprintf(fp1,'%22.14e\n',real(test_fu2(n)));
    fprintf(fp1,'%22.14e\n',imag(test_fu2(n)));
end
fclose(fp1);
fp1=fopen('test_fu3.dat','wt');
for n = 1 : changdu
    fprintf(fp1,'%22.14e\n',real(test_fu3(n)));
    fprintf(fp1,'%22.14e\n',imag(test_fu3(n)));
end
fclose(fp1);



return;
t_s = -tau_s/2 : 1/fs : tau_s/2 - 1/fs;
k_s = b/tau_s;

S=exp(1j*(pi*k_s*t_s.^2)); %���Ե�Ƶ�ź�
x_ppxs=fft(S,fft_changdu);%ƥ��ϵ��




