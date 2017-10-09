clear all;
% close all;
clc;

K=B/T;                              % 调制斜率
N = 5000;
num = 10;
% num = 10;


t=linspace(-T/2,T/2,N);
freq=linspace(-Fs/2,Fs/2,N);
St=exp(1i*pi*K*t.^2);               % 线性调频信号
coe = fft(St,N);                    % 脉压系数

fft_St = fft(St,N);                 % 纯信号脉压数据
mul_fft_St = conj(fft_St).*coe;
pc_St = fft(mul_fft_St);
pc_St_shift = fftshift(pc_St);
%%%%%%%%%%%%%%%%%%%%%%%射频噪声干扰产生%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaA = 30;                    % 干扰相对阵面的入射角度/度（从右至左）

x_sinc = -45/pi:0.05/pi:45/pi;
y_sinc = abs(sinc(x_sinc));
y_sinc_db = 20*log10(abs(y_sinc)/max(abs(y_sinc)));

xishu_A = y_sinc_db(thetaA*10+1)+13.26;   % 将第一副瓣归一化为0

Ab = 1e6;                   % 辅助通道进来的噪声干扰幅度
Aa = 10^(xishu_A/10)*Ab;    % 主通道进来的噪声干扰幅度

thetaHA = thetaA*pi/180;
juli1=zhenmian_chang*fuzhutianxianbili;
juli2=zhenmian_chang*(1-fuzhutianxianbili);
xiangweic1A=cos(thetaHA)*juli1;
xiangweic2A=cos(thetaHA)*juli2;
time1A=xiangweic1A/lamda;
time2A=xiangweic2A/lamda;
tt=1/ff;
cha1A=(time1A/tt)-fix(time1A/tt);
cha2A=(time2A/tt)-fix(time2A/tt);
xiangweicha1A=pi*cha1A;
xiangweicha2A=pi*cha2A;

% pA=fft(fir1(N-1,detlf/Fs*2));                   % 滤波器频谱
% xnA=ifft(fft(random('Normal',0,1,1,N)).*pA);    % 高斯白噪声通过滤波器
% xnA = 0.5e3 * xnA;                              % 放大噪声波形
% sum(1,N)=0;
% for i=1:N-1;
%     sum(i+1)=xnA(i)+sum(i);
% end
% xnA=sum/N;
% fp=fopen('xnA.dat','wt');
% fprintf(fp,'%s\n',1e-5);
% for i=2:length(xnA)
%     fprintf(fp,'%s\n',xnA(i));
% end;
% fclose(fp);
% return;

load xnA.dat;
xnA = xnA.';

NoiseA=Aa*(1+m*(cos(2*pi*xnA)+1j*sin(2*pi*xnA))).*cos(2*pi*f1*t);
NoiseA1=Ab*(1+m*(cos(2*pi*xnA)+1j*sin(2*pi*xnA))).*cos(2*pi*f1*t-xiangweicha1A);
NoiseA2=Ab*(1+m*(cos(2*pi*xnA)+1j*sin(2*pi*xnA))).*cos(2*pi*f1*t);
NoiseA3=Ab*(1+m*(cos(2*pi*xnA)+1j*sin(2*pi*xnA))).*cos(2*pi*f1*t+xiangweicha2A);

% figure;
% plot(real(xnA));
% title('射频噪声未调制前时域波形');
% plotbrowser('toggle');
% % return;

% figure, plot(real(NoiseA),'b');
% hold on, plot(real(NoiseA1),'r');
% hold on, plot(real(NoiseA2),'m');
% hold on, plot(real(NoiseA3),'c');hold off;
% plotbrowser('toggle');
% return;
%%%%%%%%%%%%%%%%%%%%%%%信号加上噪声干扰%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_Noise = St*1 + NoiseA*0;  % 模拟仿真数据                 % 带有信号的主通道数据
S_Noise1 = St*1 + NoiseA1*0;              % 辅助通道1数据

fft_sig = fft(S_Noise,N);
mul_fft = conj(fft_sig).*coe;
pc_result = fft(mul_fft);
pc_result_shift = fftshift(pc_result);
pc_result_shift = pc_result_shift/1e10;

fft_sig1 = fft(S_Noise1,N);
mul_fft1 = conj(fft_sig1).*coe;
pc_result1 = fft(mul_fft1);
pc_result_shift1 = fftshift(pc_result1);
pc_result_shift1 = pc_result_shift1/1e10;
%%%%%%%%%%%%%%%%%%%%%%旁瓣对消%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load 'SLC_104j3.mat'; % 实测数据
% yz = ch1(1,:)./1e6;                           % 主通道工作数据
% xg1 = slc1(1,:).*slc1_coef./1e6;              % 辅助通道1工作数据


yz = pc_result_shift;                           % 主通道工作数据
xg1 = pc_result_shift1;                         % 辅助通道1工作数据

x11 = x1';
x22 = x2';
x33 = x3';

A = zeros(3,length(x1));
A(1,:) = x1;
A(2,:) = x2;
A(3,:) = x3;

B = zeros(length(x11),3);
B(:,1) = x11;
B(:,2) = x22;
B(:,3) = x33;

C = (A*y');% 超哥互相关方法
C = C./1e14;


Ag = zeros(3,length(xg1));
Ag(1,:) = xg1;

fuzhu_zxg = A * B;
% fuzhu_zxg1 = A * A';
fuzhu_zxg = fuzhu_zxg./1e14;
% fuzhu = rand(3,num)+rand(3,num)*1j;
% fuzhu_zxg = fuzhu*fuzhu';

% fuzhu_zxg_shuchu=zeros(18,1);
% fuzhu_zxg_shuchu(1:2:end)=real(fuzhu_zxg1);
% fuzhu_zxg_shuchu(2:2:end)=imag(fuzhu_zxg1);
% fp1=fopen('fuzhu_zxg.dat','wt');          
% fwrite(fp1, fuzhu_zxg_shuchu,'float32' );
% fclose(fp1);

% zhi = det(fuzhu_zxg);
% fuzhu_zxg_ni = inv(fuzhu_zxg);
% W = fuzhu_zxg_ni * C;
W = fuzhu_zxg \ C;  
E = W' * Ag;
R = yz - E;

figure,
subplot(2,1,1),plot((abs(R))),grid,title('对消后主通道数据');
subplot(2,1,2),plot((abs(yz))),grid,title('对消前主通道数据');
% figure,plot(real(yz));
% hold on,plot(real(xg1),'r');
% plotbrowser('toggle');
% figure,plot(imag(yz));
% hold on,plot(imag(xg1),'r');
% plotbrowser('toggle');
% figure,plot(real(E));
% hold on,plot(real(R),'r');
% hold on,plot(real(pc_St_shift),'g');
% title('旁瓣对消结果实部比较');
% legend('干扰特性','对消结果','信号脉压结果');
% plotbrowser('toggle');
% figure,plot(imag(E));
% hold on,plot(imag(R),'r');
% hold on,plot(imag(pc_St_shift),'g');
% title('旁瓣对消结果虚部比较');
% legend('干扰特性','对消结果','信号脉压结果');
% plotbrowser('toggle');
%%%%%%%%%%%%%%%%%%%%%%生成测试数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fp1=fopen('_zhu_in.dat','wt');          % 主通道工作数据：   _zhu_in
% for i=1:length(yz)
%     fprintf(fp1,'%s,\n',real(pc_result_shift(i)));
%     fprintf(fp1,'%s,\n',imag(pc_result_shift(i)));
% end;
% fclose(fp1);
% fp2=fopen('_fuzhu1_in.dat','wt');       % 辅助通道1工作数据：_fuzhu1_in
% for i=1:length(xg1)
%     fprintf(fp2,'%s,\n',real(pc_result_shift1(i)));
%     fprintf(fp2,'%s,\n',imag(pc_result_shift1(i)));
% end;
%%%%%%%%%%%%%%%%%%%%%%矩阵求逆测试数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abc1 = [1+1i 2+1i 3+1i;3+1i 2+1i 1+1i;1+1i 3+1i 2+1i];
% abc = [1 1 2 1 3 1 3 1 2 1 1 1 1 1 3 1 2 1];
% det_abc = det(abc1);
% det_abc_ds = 1/det_abc;
% abc_bansui = zeros(3,3);
% abc_bansui(1,1) = det([2+1i 1+1i;3+1i 2+1i]);
% abc_bansui(1,2) = -det([2+1i 3+1i;3+1i 2+1i]);
% abc_bansui(1,3) = det([2+1i 3+1i;2+1i 1+1i]);
% abc_bansui(2,1) = -det([3+1i 1+1i;1+1i 2+1i]);
% abc_bansui(2,2) = det([1+1i 3+1i;1+1i 2+1i]);
% abc_bansui(2,3) = -det([1+1i 3+1i;3+1i 1+1i]);
% abc_bansui(3,1) = det([3+1i 2+1i;1+1i 3+1i]);
% abc_bansui(3,2) = -det([1+1i 2+1i;1+1i 3+1i]);
% abc_bansui(3,3) = det([1+1i 2+1i;3+1i 2+1i]);
% abc_ni1 = abc_bansui.*det_abc_ds;
% abc_ni2 = inv(abc1);
% 
% fp5=fopen('_fuzhu_zxg.dat','wt');       % 矩阵求逆测试数据:_fuzhu_zxg
% for i=1:18
%     fprintf(fp5,'%e\n',abc(i));
% end;
% fclose(fp5);

% dsp_fuzhu_zxg1(1,:) = '480C6909';
% dsp_fuzhu_zxg1(2,:) = '00000000';
% dsp_fuzhu_zxg1(3,:) = '47C68250';
% dsp_fuzhu_zxg1(4,:) = 'C682A0CC';
% dsp_fuzhu_zxg1(5,:) = 'C805D916';
% dsp_fuzhu_zxg1(6,:) = '460241AE';
% dsp_fuzhu_zxg1(7,:) = '47C68250';
% dsp_fuzhu_zxg1(8,:) = '4682A0CC';
% dsp_fuzhu_zxg1(9,:) = '48216E0E';
% dsp_fuzhu_zxg1(10,:) = '00000000';
% dsp_fuzhu_zxg1(11,:) = 'C80C1EBC';
% dsp_fuzhu_zxg1(12,:) = 'C61CF850';
% dsp_fuzhu_zxg1(13,:) = 'C805D916';
% dsp_fuzhu_zxg1(14,:) = 'C60241AE';
% dsp_fuzhu_zxg1(15,:) = 'C80C1EBC';
% dsp_fuzhu_zxg1(16,:) = '461CF850';
% dsp_fuzhu_zxg1(17,:) = '481647AA';
% dsp_fuzhu_zxg1(18,:) = '00000000';
% dsp_fuzhu_zxg = zeros(18,1);
% for j = 1:18
%     w1 = dsp_fuzhu_zxg1(j,:);
%     w2 = hex2dec(w1);
%     fuhaowei = bitget(uint64(w2),32);
%     if fuhaowei == 0
%         fuhao = 1;
%     elseif fuhaowei ==1
%         fuhao = -1;
%     end
%     xxx1 = bitget(w2,31:-1:24);
%     yyy1 = bitget(w2,23:-1:1);
%     xxx2 = 0;
%     for i = 1:8
%         xxx3 = xxx1(1,i)*(2^(8-i));
%         xxx2 = xxx3 + xxx2;
%     end
%     yyy2 = 0;
%     for i = 1:23
%         yyy3 = yyy1(1,i)*(2^(23-i));
%         yyy2 = yyy3 + yyy2;
%     end
%     dsp_fuzhu_zxg(j,1) = fuhao * (1 + yyy2/2^23)*2^(xxx2-127);
% end
% aaa = zeros(3,3);
% aaa(1,1) = dsp_fuzhu_zxg(1,1) + dsp_fuzhu_zxg(2,1)*1i;
% aaa(1,2) = dsp_fuzhu_zxg(3,1) + dsp_fuzhu_zxg(4,1)*1i;
% aaa(1,3) = dsp_fuzhu_zxg(5,1) + dsp_fuzhu_zxg(6,1)*1i;
% aaa(2,1) = dsp_fuzhu_zxg(7,1) + dsp_fuzhu_zxg(8,1)*1i;
% aaa(2,2) = dsp_fuzhu_zxg(9,1) + dsp_fuzhu_zxg(10,1)*1i;
% aaa(2,3) = dsp_fuzhu_zxg(11,1) + dsp_fuzhu_zxg(12,1)*1i;
% aaa(3,1) = dsp_fuzhu_zxg(13,1) + dsp_fuzhu_zxg(14,1)*1i;
% aaa(3,2) = dsp_fuzhu_zxg(15,1) + dsp_fuzhu_zxg(16,1)*1i;
% aaa(3,3) = dsp_fuzhu_zxg(17,1) + dsp_fuzhu_zxg(18,1)*1i;

% aaa = fuzhu_zxg;
% 
% bbb = det(aaa);
% ccc1 = aaa(1,1)*aaa(2,2)*aaa(3,3)+aaa(1,3)*aaa(2,1)*aaa(3,2)+aaa(1,2)*aaa(2,3)*aaa(3,1);
% ccc2 = -aaa(1,3)*aaa(2,2)*aaa(3,1)-aaa(1,2)*aaa(2,1)*aaa(3,3)-aaa(1,1)*aaa(2,3)*aaa(3,2);
% ccc = ccc1 + ccc2;
% 
% % bbb = -268435456 - 16777216*1i;
% zhi_daoshu = (real(ccc)-imag(ccc)*1i)/(real(ccc)^2+imag(ccc)^2);
% 
% % v1 = real(bbb)^2;
% % v2 = imag(bbb)^2; 
% % v3 = real(bbb)^2+imag(bbb)^2;
% % v4 = 1/(real(bbb)^2+imag(bbb)^2);
% % v5 = real(bbb)*(1/(real(bbb)^2+imag(bbb)^2));
% % v6 = -imag(bbb)*(1/(real(bbb)^2+imag(bbb)^2));
% 
% fuzhu_zxg_bs = zeros(3,3);
% fuzhu_zxg_bs(1,1) = aaa(2,2)*aaa(3,3)-aaa(2,3)*aaa(3,2);
% fuzhu_zxg_bs(1,2) = -(aaa(1,2)*aaa(3,3)-aaa(1,3)*aaa(3,2)); 
% fuzhu_zxg_bs(1,3) = aaa(1,2)*aaa(2,3)-aaa(1,3)*aaa(2,2);
% fuzhu_zxg_bs(2,1) = -(aaa(2,1)*aaa(3,3)-aaa(2,3)*aaa(3,1));
% fuzhu_zxg_bs(2,2) = aaa(1,1)*aaa(3,3)-aaa(1,3)*aaa(3,1);
% fuzhu_zxg_bs(2,3) = -(aaa(1,1)*aaa(2,3)-aaa(1,3)*aaa(2,1));
% fuzhu_zxg_bs(3,1) = aaa(2,1)*aaa(3,2)-aaa(2,2)*aaa(3,1);
% fuzhu_zxg_bs(3,2) = -(aaa(1,1)*aaa(3,2)-aaa(1,2)*aaa(3,1));
% fuzhu_zxg_bs(3,3) = aaa(1,1)*aaa(2,2)-aaa(1,2)*aaa(2,1);
% 
% zhi_daoshu = -3.718014340847731e-09 - 2.323758963029832e-10i;
% fuzhu_zxg_ni1 = fuzhu_zxg_bs.*zhi_daoshu;

% aaa = fuzhu_zxg;% 矩阵除法的另一种方法：高斯消去法
% bbb = C;
% aaa1(1,:) = aaa(1,:);
% bbb1(1,:) = bbb(1,:);
% aaa1(2,:) = aaa(2,:)-aaa(1,:).*(aaa(2,1)/aaa(1,1));
% bbb1(2,:) = bbb(2,:)-bbb(1,:).*(aaa(2,1)/aaa(1,1));
% aaa1(3,:) = aaa(3,:)-aaa(1,:).*(aaa(3,1)/aaa(1,1));
% bbb1(3,:) = bbb(3,:)-bbb(1,:).*(aaa(3,1)/aaa(1,1));
% aaa2(1,:) = aaa1(1,:);
% bbb2(1,:) = bbb1(1,:);
% aaa2(2,:) = aaa1(2,:);
% bbb2(2,:) = bbb1(2,:);
% aaa2(3,:) = aaa1(3,:)-aaa1(2,:).*(aaa1(3,2)/aaa1(2,2));
% bbb2(3,:) = bbb1(3,:)-bbb1(2,:).*(aaa1(3,2)/aaa1(2,2));
% www(3,1) = bbb2(3,1)/aaa2(3,3);
% www(2,1) = (bbb2(2,1)-www(3,1)*aaa2(2,3))/aaa2(2,2);
% www(1,1) = (bbb2(1,1)-www(3,1)*aaa2(1,3)-www(2,1)*aaa2(1,2))/aaa2(1,1);

ccc = fuzhu_zxg;% 矩阵除法的另一种方法：高斯消去法
ddd = C;
ccc1(1,:) = ccc(1,:);
ddd1(1,:) = ddd(1,:);
ccc1(2,:) = ccc(2,:).*ccc(1,1)-ccc(1,:).*ccc(2,1);
ddd1(2,:) = ddd(2,1).*ccc(1,1)-ddd(1,1).*ccc(2,1);
ccc1(3,:) = ccc(3,:).*ccc(1,1)-ccc(1,:).*ccc(3,1);
ddd1(3,:) = ddd(3,1).*ccc(1,1)-ddd(1,1).*ccc(3,1);
ccc2(1,:) = ccc1(1,:);
ddd2(1,:) = ddd1(1,:);
ccc2(2,:) = ccc1(2,:);
ddd2(2,:) = ddd1(2,:);
ccc2(3,:) = ccc1(3,:).*ccc1(2,2)-ccc1(2,:).*ccc1(3,2);
ddd2(3,:) = ddd1(3,:).*ccc1(2,2)-ddd1(2,:).*ccc1(3,2);
wwww(3,1) = ddd2(3,1)/ccc2(3,3);
wwww(2,1) = (ddd2(2,1)-wwww(3,1)*ccc2(2,3))/ccc2(2,2);
wwww(1,1) = (ddd2(1,1)-wwww(3,1)*ccc2(1,3)-wwww(2,1)*ccc2(1,2))/ccc2(1,1);

a = real(fuzhu_zxg(1,1));
c = real(fuzhu_zxg(1,2));
d = imag(fuzhu_zxg(1,2));
e = real(fuzhu_zxg(1,3));
f = imag(fuzhu_zxg(1,3));
g = real(fuzhu_zxg(2,2));
h = imag(fuzhu_zxg(2,2));
m = real(fuzhu_zxg(2,3));
n = imag(fuzhu_zxg(2,3));
o = real(fuzhu_zxg(3,3));

% x4=a*o-e^2-f^2;
% x3=a*n-c*f+d*e;
% x2=a*m-c*e-d*f;
% x1=g*a-c^2-d^2;
% x5=x1*x4-x2^2-x3^2;
% x6 = a*(o*(a*g-c*c-d*d)-g*(e*e+f*f)-a*(m*m+n*n)+2*c*(e*m+f*n)+2*d*(f*m-e*n));
x6 = (o*(a*g-c*c-d*d)-g*(e*e+f*f)-a*(m*m+n*n)+2*c*(e*m+f*n)+2*d*(f*m-e*n));

x1 = real(C(1,1));
y1 = imag(C(1,1));
x2 = real(C(2,1));
y2 = imag(C(2,1));
x3 = real(C(3,1));
y3 = imag(C(3,1));

z2 = (a*x2-c*x1-d*y1)+(a*y2-c*y1+d*x1)*1j;
zz32 = (a*m-c*e-d*f)+(c*f-a*n-d*e)*1j;
zz22 = a*g-c^2-d^2;
z31 = (a*x3-e*x1-f*y1)+(a*y3-e*y1+f*x1)*1j;
z32 = z31*zz22-z2*zz32;

zzz3real = ((a*g-c^2-d^2)*x3+(c*m-d*n-e*g)*x1+(d*m+c*n-f*g)*y1+(c*e+d*f-a*m)*x2+(c*f-a*n-d*e)*y2);
zzz3imag = ((a*g-c^2-d^2)*y3+(f*g-d*m-c*n)*x1+(c*m-e*g-d*n)*y1+(a*n+d*e-c*f)*x2+(c*e+d*f-a*m)*y2);





