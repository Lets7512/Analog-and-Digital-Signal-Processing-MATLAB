bit_rate = 64;
bit_width = 100;
T=2;
data_polar = [];
data_Manchester = [];
Voltage = 1;
for i=1:T*bit_rate
    r = randi([0 1],1,1);
    r_d = 1*(r==1) + -1*(r==0);
    d_p = r_d*ones(1,bit_width);
    data_polar = [data_polar , d_p];
    d_low_high = [-1*ones(1,bit_width/2),ones(1,bit_width/2)];
    d_high_low = [ones(1,bit_width/2),-ones(1,bit_width/2)];
    d_m = d_high_low*(r==1) + d_low_high*(r==0);
    data_Manchester = [data_Manchester , d_m];
end
data_polar=data_polar*Voltage;
data_Manchester=data_Manchester*Voltage;
t = linspace(0,T,length(data_polar));
figure;
plot(t,data_polar,'LineWidth',1);
title('Line code Polar');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
figure;
plot(t,data_Manchester,'LineWidth',1);
title('Line code MANCHESTER');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
%-----------------------------------------------------------
%-----------------------------------------------------------
delta_t = (t(1,2)-t(1,1));
Fs = bit_rate*bit_width;
L = length(t);
figure_accuracy = 20;
nfft = 2^nextpow2(L)*figure_accuracy;
f = linspace(-Fs/2,Fs/2,nfft);
f_freq = linspace(-bit_rate*2,2*bit_rate,nfft);
bit_time = 1/bit_rate;
N = (bit_time/delta_t);
figure;
FT_Polar = fft(data_polar,nfft);
FT_Polar = fftshift(FT_Polar);
FT_Polar = abs(FT_Polar)/Fs;
plot(f_freq,FT_Polar);
title('Fast Fourier Transform of Polar line coding');
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
legend("Polar");
%-----------------------------------------------------------
%-----------------------------------------------------------
figure;
FT_Manchester = fft(data_Manchester,nfft);
FT_Manchester = fftshift(FT_Manchester);
FT_Manchester = abs(FT_Manchester)/Fs;
plot(f_freq,FT_Manchester);
title('Fast Fourier Transform of Manchester line coding');
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
legend("Manchester");
%-----------------------------------------------------------
%-----------------------------------------------------------
one_bit_polar_data = data_polar(1,1:bit_width);
FT_one_bit_polat_data = fft(one_bit_polar_data,nfft);
FT_one_bit_polat_data = fftshift(FT_one_bit_polat_data);
FT_one_bit_polat_data = abs(FT_one_bit_polat_data);
figure;
p_p = FT_one_bit_polat_data;
p_p = p_p.^2/Fs/N;
plot(f(f>=0 & f<=2*bit_rate),p_p(f>=0 & f<=2*bit_rate));
title('Power Spectrum of Polar line coding');
xlabel('Frequency (Hz)');
ylabel('V(f)^2');
legend("Polar");
%-----------------------------------------------------------
%-----------------------------------------------------------
one_bit_man_data = data_Manchester(1,1:bit_width);
FT_one_bit_man_data = fft(one_bit_man_data,nfft);
FT_one_bit_man_data = fftshift(FT_one_bit_man_data);
FT_one_bit_man_data = abs(FT_one_bit_man_data);
figure;
man = FT_one_bit_man_data;
man = man.^2/Fs/N;
plot(f(f>=0 & f<=2*2*bit_rate),man(f>=0 & f<=2*2*bit_rate));
title('Power Spectrum of Manchester line coding');
xlabel('Frequency (Hz)');
ylabel('V(f)^2');
legend("Manchester");
%-----------------------------------------------------------
%-----------------------------------------------------------
figure;
plot(f(f>=0 & f<=2*2*bit_rate),p_p(f>=0 & f<=2*2*bit_rate));
hold on
plot(f(f>=0 & f<=2*2*bit_rate),man(f>=0 & f<=2*2*bit_rate));
title(sprintf('Power Spectrum of Manchester vs Polar where R{b} = %d , V = %d',bit_rate,Voltage));
xlabel('Frequency (Hz)');
ylabel('V(f)^2');
legend("Manchester");
legend("Polar","Manchester")