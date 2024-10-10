clc; clear; close all; 
  
sigma_bp = 1/2; 
 
PtDbm = 10;   
Pt = 10.^(PtDbm/10)*10^-3; % Công suất phát 
 
NDbm = -20;  
N  = 10.^(NDbm/10)*10^-3; % Công suất nhiễu AWGN 
  
PminDbm = -20:-5;   
Pmin = 10.^(PminDbm/10)*10^-3; % Ngưỡng dịch vụ 
  
d_ch = 2;      % Khoảng cách chuẩn hóa 
alpha = 2.5;   % Hệ số suy hao đường truyền 
 
% ================ Mô phỏng ================= 
Ntry = 10^6; 
OP_simul = zeros(1,length(Pmin)); 
for i = 1:length(OP_simul) 
  
    % ========= Phần thực, phần ảo ========== 
    hI = sqrt(sigma_bp) * randn(1,Ntry); 
    hQ = sqrt(sigma_bp) * randn(1,Ntry); 
  
    % ======= Công suất thu tức thời ======== 
    Pr = Pt*d_ch.^-alpha.*(hI.^2 + hQ.^2) + N;          % Công suất thu tức thời 
    indexLess = find(Pr < Pmin(i));   % Tìm vị trí dừng 
     
    OP_simul(i) = length(indexLess)/Ntry; % Xác suất dừng mô phỏng 
end 
OP_simul; 
  
% ============== Lý thuyết ================== 
a = (Pmin - N)/(Pt*d_ch.^-alpha);
OP_theo = 1- exp(-a/(2*sigma_bp));    % Xác suất dừng lý thuyết 
  
figure(1) 
plot(PminDbm,OP_simul,'ko','linewidth',1.4);  
hold on; 
plot(PminDbm,OP_theo,'r-','linewidth',1.4); 
xlabel('Pmin (dBm)'); ylabel('OP'); 
legend('Simulation', 'Theory');