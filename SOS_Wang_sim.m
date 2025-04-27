% SOS_Wang_sim: This program simulate the statistical properties of SOS
% -based Rayleigh fading chaneel using Wang Model, 
% DOI: https://doi.org/10.1007/s12209-012-1888-1
% Copyright (C) 2025  Mohammad Safa
% GitHub Repository: https://github.com/mhr98/SOS-Wang-sim
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%% generate rayleigh parameters
fd = 100; % maximum Doppler shift
Ts = 1e-4; % smple duration

L = 64000; % number of samples for simulation

M= 2^12;  %fft samples
t= [0:M-1]*Ts;
f= [-M/2:M/2-1]/(M*Ts*fd);

N = 8; % number of sinusiod 

%% generationg random values

theta=(rand-0.5)*2*pi;
for n=1:N

    phi_i(n) = (rand-0.5);
    alpha_i(n) = fd * Ts * sin((2*pi*n + theta/32 - pi*(3/2))/(4*N));
    
    phi_j(n) = (rand-0.5);
    alpha_j(n) = fd * Ts * cos((2*pi*n + theta/32 - pi*(3/2))/(4*N));
end

% printing (optional)

% fprintf('phi_i:\n');
% for j = 1:N
%     
%         fprintf('%.7f', phi_i(j));
%         fprintf(',');
% end
% fprintf('\n\n');
% 
% fprintf('alpha_i:\n');
% for j = 1:N
%     
%         fprintf('%.7f', alpha_i(j));
%         fprintf(',');
% end
% fprintf('\n\n');
% 
% fprintf('phi_j:\n');
% for j = 1:N
%     
%         fprintf('%.7f', phi_j(j));
%         fprintf(',');
% end
% fprintf('\n\n');
% 
% fprintf('alpha_j:\n');
% for j = 1:N
%     
%         fprintf('%.7f', alpha_j(j));
%         fprintf(',');
% end
% fprintf('\n\n');


%% sos calaculations

r=[];

alpha_i_acc = zeros(1,N);
alpha_j_acc = zeros(1,N);

for i = 1: L
    
    si = 0; 
    sj = 0;

    for n=1:N
        
        alpha_i_acc(n) = alpha_i_acc(n) + alpha_i(n);
        si = si + sin(2*pi*(alpha_i_acc(n)+phi_i(n)));
        
        alpha_j_acc(n) = alpha_j_acc(n) + alpha_j(n);
        sj = sj + sin(2*pi*(alpha_j_acc(n)+phi_j(n)));
    end
    
    r = [r sqrt(2/N)*(si + 1i*sj)];
    
end

%% Autocorrelation Calculatio

Ns = L;
temp=zeros(2,Ns);
for i=1:Ns
   j=i:Ns; 
   temp(1:2,j-i+1)= temp(1:2,j-i+1)+[r(i)'*r(j); ones(1,Ns-i+1)];
end
k=1:M; 
Simulated_corr(k)= temp(1,k)./temp(2,k);
Classical_corr= 2*besselj(0,2*pi*fd*t);  %theoritical

% Fourier transform of autocorrelation
Classical_Y= fftshift(fft(Classical_corr));
Simulated_Y= fftshift(fft(Simulated_corr));

%% plotting the results

figure(1);
hold on;grid on; box on;
xlabel('time[s]','FontSize',18,'FontName', 'Times'); 
ylabel('Magnitude[dB]','FontSize',18,'FontName', 'Times'); 
plot([1:L]*Ts,10*log10(abs(r)),'LineWidth',2); 
axis([0 (L/10)*Ts 1.1*min(10*log10(abs(r))) 1.1*max(10*log10(abs(r)))])
set(gca,'fontsize',18);

figure(2);
hold on;grid on; box on;
xlabel('Magnitude','FontSize',18,'FontName', 'Times'); 
ylabel('PDF','FontSize',18,'FontName', 'Times');
[pdfx, xval]=hist(abs(r),30);
bar(xval,pdfx/(sum(pdfx)*(xval(2)-xval(1))));
set(gca,'fontsize',18);

figure(3);
hold on;grid on; box on;
xlabel('Phase[rad]','FontSize',18,'FontName', 'Times'); 
ylabel('PDF','FontSize',18,'FontName', 'Times');
[pdfx, xval]=hist(angle(r),30);
bar(xval,pdfx/(sum(pdfx)*(xval(2)-xval(1))));
set(gca,'fontsize',18);

figure(4);
hold on;grid on; box on;
xlabel('delay \tau [s]','FontSize',18,'FontName', 'Times'); 
ylabel('Correlation','FontSize',18,'FontName', 'Times');
plot(t,abs(Classical_corr),'b-','LineWidth',2.5)
plot(t,abs(Simulated_corr),'r:','LineWidth',3.5);
legend({'Theoretical','Simulated'},'FontSize',18,'FontName', 'Times')
xlim([0 0.2]);
set(gca,'fontsize',18);

figure(5);
hold on;grid on; box on;
xlabel('f/f_d','FontSize',18,'FontName', 'Times'); 
ylabel('Magnitude','FontSize',18,'FontName', 'Times'); 
plot(f,abs(Classical_Y),'b-','LineWidth',2.5)
plot( f,abs(Simulated_Y),'r:','LineWidth',3.5);
axis([-1 1 0 1.1*max(abs(Classical_Y))])
legend({'Theoretical','Simulated'},'FontSize',18,'FontName', 'Times')
set(gca,'fontsize',18);



