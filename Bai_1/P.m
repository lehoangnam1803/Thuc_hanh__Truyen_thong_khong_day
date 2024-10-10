clc; clear; close all;

f = 900*10^6; % Tần số phát
c = 3*10^8;
lamda = c/f;  % Bước sóng phát
d0 = 1;       % Không các tham chiếu

K = 20*log10(lamda/(4*pi*d0));     % Hệ số K

% ========= Dữ liệu đo đạt =============
d = [5 25 65 110 400 1000];
G = [-60 -80 -105 -115 -135 -150];
N = length(d);

% =======Khớp phương trình bằng MMSE ===
syms alpha;
Falpha = 0;
for i = 1:N
    Hi = K - 10*alpha*log10(d(i)/d0);
    Falpha = Falpha + (G(i) - Hi).^2;
end

Falpha = 1/N * Falpha;  % Hàm MMSE
Falpha = simplify(Falpha);

% ========= Tìm alpha tối ưu ===========
Falpha_diff = diff(Falpha,alpha);

alpha_opt = vpasolve(Falpha_diff == 0,alpha);
alpha_opt = double(alpha_opt)

% == Phương sai cho mô hình Shadowing ==
sigmaPsiDb = subs(Falpha, alpha_opt); % thay giá trị alpha tìm được vào hàm f(alpha)
sigmaPsiDb = double(sigmaPsiDb)

%% Sử dụng mô hình xây dựng ở trên
% thực hiện dự đoán suy hao tại các khoảng cách

d_dudoan = [50 150 200 300];
Hi_dudoan = K - 10*alpha_opt*log10(d_dudoan./d0)