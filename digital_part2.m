%------------------------------------------------------------
%------------------------------------------------------------
%                     Digital Transmitter
%------------------------------------------------------------
%------------------------------------------------------------
%---------------------------------------------------------
%---------------------------------------------------------
bit_rate = 4;
T=2;
bit_width = 100;
%resolution = 1/(bit_rate*figure_accuracy);
A=3;
F1=bit_rate*4;
F3=bit_rate*4;
F2=bit_rate*2;
%------------------------------------------------------------
%------------------------------------------------------------
data_ask = [];
data_fsk = [];
data_psk = [];

for i=1:T*bit_rate
    r = randi([0 1],1,1);
    r_ask = r;
    r_fsk = r;
    r_psk = r*(r==1) + -1*(r==0);
    
    d_ask = r_ask*ones(1,bit_width);
    d_fsk = r_fsk*ones(1,bit_width);
    d_psk = r_psk*ones(1,bit_width);
    
    data_ask = [data_ask,d_ask];
    data_fsk = [data_fsk,d_fsk];
    data_psk = [data_psk,d_psk];
end
t = linspace(0,T,length(data_psk));
figure;
subplot(3,1,1);
plot(t,data_ask,'LineWidth',2)
title('ASK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Binary Data");
subplot(3,1,2);
plot(t,data_fsk,'LineWidth',2)
title('FSK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Binary Data");
subplot(3,1,3);
plot(t,data_psk,'LineWidth',2)
title('PSK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Binary Data");
%------------------------------------------------------------
%------------------------------------------------------------
%ASK
Cos=A.*cos(2*pi*F1*t);
Sin=A.*sin(2*pi*F1*t);
ASK=Cos;
ASK_d = @(b) ASK.*(b==1);
figure;
subplot(2,1,1);plot(t,ASK_d(data_ask))
title('Data as ASK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("ASK Data");
subplot(2,1,2);plot(t,data_ask,'LineWidth',2)
title('ASK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Binary Data");
%------------------------------------------------------------
%------------------------------------------------------------
%FSK
sin_alph=A.*cos(2*pi*F3*t);
sin_beta=A.*cos(2*pi*F2*t);
%FSK=sin_alph.*(A/2.*square(2*pi*F2*t+pi)+A/2)+ASK;
FSK_d = @(b) sin_alph.*(b==1) + sin_beta.*(b==0);
figure;
subplot(2,1,1);plot(t,FSK_d(data_fsk))
title('Data as FSK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("ASK Data");
subplot(2,1,2);plot(t,data_fsk,'LineWidth',2)
title('FSK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Binary Data");
%------------------------------------------------------------
%------------------------------------------------------------
%PSK
alpha=pi;
PSK=A.*cos(2*pi*F1*t+alpha);
%PSK=ASK;
PSK_d = @(b) Cos.*(b==1) + PSK.*(b==-1);
figure;
subplot(2,1,1);plot(t,PSK_d(data_psk))
title('Data as PSK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("PSK Data");
subplot(2,1,2);plot(t,data_psk,'LineWidth',2)
title('PSK');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Binary Data");
%------------------------------------------------------------
%------------------------------------------------------------
figure;
Fs = bit_rate*bit_width;
L = length(t);
nfft = 2^nextpow2(L)*bit_width;
f = linspace(-Fs/2,Fs/2,nfft);
ask_data = ASK_d(data_ask);
spec_ask = fft(ask_data,nfft);
spec_ask = fftshift(spec_ask);
spec_ask = abs(spec_ask);
spec_ask = spec_ask/Fs;
plot(f,spec_ask);
title('Transmitted ASK Spectrum Response');
ylabel('Voltage (V)');
xlabel('Frequency (Hz)');
legend("ASK");
%------------------------------------------------------------
%------------------------------------------------------------
figure;
fsk_data = FSK_d(data_fsk);
spec_fsk = fft(fsk_data,nfft);
spec_fsk = abs(spec_fsk);
spec_fsk = fftshift(spec_fsk);
spec_fsk = spec_fsk/Fs;
plot(f,spec_fsk);
title('Transmitted FSK Spectrum Response');
ylabel('Voltage (V)');
xlabel('Frequency (Hz)');
legend("FSK");
%------------------------------------------------------------
%------------------------------------------------------------
figure;
psk_data = PSK_d(data_psk);
spec_psk = fft(psk_data,nfft);
spec_psk = abs(spec_psk);
spec_psk = fftshift(spec_psk);
spec_psk = spec_psk/Fs;
plot(f,spec_psk);
title('Transmitted PSK Spectrum Response');
ylabel('Voltage (V)');
xlabel('Frequency (Hz)');
legend("PSK");
%------------------------------------------------------------
%------------------------------------------------------------
%                     Digital Reciever
%------------------------------------------------------------
%------------------------------------------------------------
%---------------------------------------------------------
%---------------------------------------------------------
%                    Before Filteration
%---------------------------------------------------------
%---------------------------------------------------------
%PSK Reciever
A=1;
phase_shift_psk=0;
coherent_signal_psk=A.*cos(2*pi*F1*t+phase_shift_psk);
sig_psk=psk_data.*coherent_signal_psk;
%reveiver_filtered = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
figure;
plot(t,psk_data,'LineWidth',1);
hold on;
plot(t,sig_psk,'LineWidth',1);
hold off;
title('Transmitted PSK Signal vs Recieved PSK Signal');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Transmitted Signal","Recieved Signal");
%------------------------------------------------------------
%------------------------------------------------------------
figure;
plot(t,sig_psk,'LineWidth',1);
title('Recieved PSK Signal Before Filter');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal");
%------------------------------------------------------------
%------------------------------------------------------------
figure;
spec_received_psk = fft(sig_psk,nfft);
spec_received_psk = abs(spec_received_psk);
spec_received_psk = fftshift(spec_received_psk);
spec_received_psk = spec_received_psk/Fs;
plot(f,spec_received_psk);
title('Recieved PSK Spectrum Response');
ylabel('Voltage (V)');
xlabel('Frequency (Hz)');
legend("PSK");
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 0
figure;
phase_shift_psk=0;
coherent_signal_psk=A.*cos(2*pi*F1*t+phase_shift_psk);
sig_psk=psk_data.*coherent_signal_psk;
plot(t,sig_psk,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
filtered_psk = lowpass(sig_psk,bit_rate,length(sig_psk));
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
plot(t,filtered_psk,'LineWidth',1);
plot(t,data_psk,'LineWidth',1);
hold off;
title(sprintf('Original PSK Data vs Recieved PSK Data Under phase effect = %d',phase_shift_psk));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 30
figure;
phase_shift_psk=30;
coherent_signal_psk=A.*cos(2*pi*F1*t+phase_shift_psk);
sig_psk=psk_data.*coherent_signal_psk;
plot(t,sig_psk,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
filtered_psk = lowpass(sig_psk,bit_rate,length(sig_psk));
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
plot(t,filtered_psk,'LineWidth',1);
plot(t,data_psk,'LineWidth',1);
hold off;
title(sprintf('Original PSK Data vs Recieved PSK Data Under phase effect = %d',phase_shift_psk));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 60
figure;
phase_shift_psk=60;
coherent_signal_psk=A.*cos(2*pi*F1*t+phase_shift_psk);
sig_psk=psk_data.*coherent_signal_psk;
plot(t,sig_psk,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
filtered_psk = lowpass(sig_psk,bit_rate,length(sig_psk));
plot(t,filtered_psk,'LineWidth',1);
plot(t,data_psk,'LineWidth',1);
hold off;
title(sprintf('Original PSK Data vs Recieved PSK Data Under phase effect = %d',phase_shift_psk));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 90
figure;
phase_shift_psk=90;
coherent_signal_psk=A.*cos(2*pi*F1*t+phase_shift_psk);
sig_psk=psk_data.*coherent_signal_psk;
plot(t,sig_psk,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
filtered_psk = lowpass(sig_psk,bit_rate,length(sig_psk));
plot(t,filtered_psk,'LineWidth',1);
plot(t,data_psk,'LineWidth',1);
hold off;
title(sprintf('Original PSK Data vs Recieved PSK Data Under phase effect = %d',phase_shift_psk));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%---------------------------------------------------------
%---------------------------------------------------------
%---------------------------------------------------------
%---------------------------------------------------------
%ASK Reciever
phase_shift_ask=0;
coherent_signal_ask=A.*cos(2*pi*F1*t+phase_shift_ask);
sig_ask=ask_data.*coherent_signal_ask;
%filter_u = @(b) 1.*(b>0) + -1.*(b<0);
%filtered = filter_u(u);
figure;
plot(t,ask_data,'LineWidth',1);
hold on;
plot(t,sig_ask,'LineWidth',1);
hold off;
title('Transmitted ASK Signal vs Recieved ASK Signal');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Transmitted Signal","Recieved Signal");
%---------------------------------------------------------
%---------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
figure;
plot(t,sig_ask,'LineWidth',1);
title('Recieved ASK Signal Before Filter');
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal");
%------------------------------------------------------------
%------------------------------------------------------------
figure;
spec_received_ask = fft(sig_ask,nfft);
spec_received_ask = abs(spec_received_ask);
spec_received_ask = fftshift(spec_received_ask);
spec_received_ask = spec_received_ask/Fs;
plot(f,spec_received_ask);
title('Recieved ASK Spectrum Response');
ylabel('Voltage (V)');
xlabel('Frequency (Hz)');
legend("ASK");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 0
figure;
phase_shift_ask=0;
coherent_signal_ask=A.*cos(2*pi*F1*t+phase_shift_ask);
sig_ask=ask_data.*coherent_signal_ask;
plot(t,sig_ask,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
filtered_ask = lowpass(sig_ask,bit_rate,length(sig_ask));
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
plot(t,filtered_ask,'LineWidth',1);
plot(t,data_ask,'LineWidth',1);
hold off;
title(sprintf('Original ASK Data vs Recieved ASK Data Under phase effect = %d',phase_shift_ask));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 30
figure;
phase_shift_ask=30;
coherent_signal_ask=A.*cos(2*pi*F1*t+phase_shift_ask);
sig_ask=ask_data.*coherent_signal_ask;
plot(t,sig_ask,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
filtered_ask = lowpass(sig_ask,bit_rate,length(sig_ask));
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
plot(t,filtered_ask,'LineWidth',1);
plot(t,data_ask,'LineWidth',1);
hold off;
title(sprintf('Original ASK Data vs Recieved ASK Data Under phase effect = %d',phase_shift_ask));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 60
figure;
phase_shift_ask=60;
coherent_signal_ask=A.*cos(2*pi*F1*t+phase_shift_ask);
sig_ask=ask_data.*coherent_signal_ask;
plot(t,sig_ask,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
filtered_ask = lowpass(sig_ask,bit_rate,length(sig_ask));
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
plot(t,filtered_ask,'LineWidth',1);
plot(t,data_ask,'LineWidth',1);
hold off;
title(sprintf('Original ASK Data vs Recieved ASK Data Under phase effect = %d',phase_shift_ask));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
%ocsilate phase = 90
figure;
phase_shift_ask=90;
coherent_signal_ask=A.*cos(2*pi*F1*t+phase_shift_ask);
sig_ask=ask_data.*coherent_signal_ask;
plot(t,sig_ask,'LineWidth',1);
hold on;
%filter_psk = @(b) 1.*(b>0) + -1.*(b<0);
%filtered_psk = filter_psk(sig_psk);
filtered_ask = lowpass(sig_ask,bit_rate,length(sig_ask));
%filtered_psk = abs((ifft(fft(sig_psk).*fft(sinc(bit_rate*2*t)))));
plot(t,filtered_ask,'LineWidth',1);
plot(t,data_ask,'LineWidth',1);
hold off;
title(sprintf('Original ASK Data vs Recieved ASK Data Under phase effect = %d',phase_shift_ask));
ylabel('Voltage (V)');
xlabel('time (sec)');
legend("Recieved Signal(before filter)","Filtered Signal","Original Data");
%------------------------------------------------------------
%------------------------------------------------------------