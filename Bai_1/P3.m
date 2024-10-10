clc; clear; close all; 
  
Ntry = 10^6; 
sigma_bp = 1/2; 
% ============ Phần thực, phần ảo =========== 
hI = sqrt(sigma_bp) * randn(1, Ntry); 
hQ = sqrt(sigma_bp) * randn(1, Ntry); 
  
% ============ Độ lớn kênh ================== 
abs_h = sqrt(hI.^2 + hQ.^2); 
  
% ============ Độ lợi kênh ================== 
abs_h_bp = hI.^2 + hQ.^2; 
  
% ========= Monte Carlo Ham Mat do PDF ====== 
[fhI,xI] = ksdensity(hI); 
[fhQ,xQ] = ksdensity(hQ); 
[fabs_h,x_abs] = ksdensity(abs_h); 
[fabs_h_bp,x_abs_bp] = ksdensity(abs_h_bp); 
  
% ========= Lý thuyết hàm mật độ PDF ======== 
fhIL = 1/sqrt(2*pi.*sigma_bp) * exp(-xI.^2/(2.*sigma_bp)); 
fhQL = 1/sqrt(2*pi.*sigma_bp) * exp(-xQ.^2/(2.*sigma_bp)); 
fabs_hL = x_abs/sigma_bp .* exp(-x_abs.^2/(2.*sigma_bp)); 
fabs_h_bpL = 1/(2.*sigma_bp) * exp(-x_abs_bp/(2.*sigma_bp)); 
  
% ============== Đồ thị ===================== 
figure(1) 
plot(xI,fhI,'ko','linewidth',1.4); hold on; 
plot(xI,fhIL,'r-','linewidth',1.4);  
xlabel('h_I'); ylabel('f(h_I)'); 
legend('Simulation','Theory') 
  
figure(2) 
plot(xQ,fhQ,'ko','linewidth',1.4); hold on; 
plot(xQ,fhQL,'b-','linewidth',1.4);  
xlabel('h_Q'); ylabel('f(h_Q)'); 
legend('Simulation','Theory') 
  
figure(3) 
plot(x_abs,fabs_h,'ko','linewidth',1.4); hold on; 
plot(x_abs,fabs_hL,'g-','linewidth',1.4);  
xlabel('|h|'); ylabel('f(|h|)'); 
legend('Simulation','Theory') 
  
figure(4) 
plot(x_abs_bp,fabs_h_bp,'ko','linewidth',1.4); 
hold on; 
plot(x_abs_bp,fabs_h_bpL,'m-','linewidth',1.4);  
xlabel('|h|^2'); ylabel('f(|h|^2)'); 
legend('Simulation','Theory')