clc; clear; close all;
f = 900*10^6;   % Tần số phát
c = 3*10^8;
lamda = c/f;    % Bước sóng phát
d0 = 1;         % Khoảng các tham chiếu

K = 20*log10(lamda/(4*pi*d0));% Hệ số K

alpha = 3.9581;
sigmaPsiDb = 9.681;

d = 100;
Pt = 10; PtDbm = 10*log10(Pt);
Pmin = -110:-100; % Ngưỡng dịch vụ

% ============= Monte Carlo ==============
Ntry = 10^5;
OP_Monte = zeros(1,length(Pmin));
for i = 1:length(Pmin)
    PsiDb = sqrt(sigmaPsiDb) * randn(1,Ntry);  % Tạo các mẫu theo Gauss
    PrDbm = PtDbm + K - 10*alpha*log10(d/d0) - PsiDb;  % Công suất thu
    indexLess = find(PrDbm < Pmin(i)); % Tìm vị trí dừng
    OP_Monte(i) = length(indexLess) / Ntry; % Xác suất dừng
end
OP_Monte

% ============== Lý thuyết ===============
a =  PtDbm + K - 10*alpha*log10(d/d0) - Pmin;
OP_Theo = qfunc(a/sqrt(sigmaPsiDb));

% =============== Đồ thị =================
figure(1)
semilogy(Pmin,OP_Theo,'r-','linewidth',1.4);
hold on;
semilogy(Pmin,OP_Monte,'ko','linewidth',1.4)
xlabel('Pmin (dBm)'); ylabel('OP');
legend('Theo','Simulation')