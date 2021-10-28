%1. Plot the function
y = @(t) (2-t).*(t>=0 & t<=2) + (2+t).*(t<0 & t>=-2) + 0.*(t>2) + 0.*(t<-2);
Fs = 10;            % Sampling frequency
resulotion = 0.01;
N=Fs/resulotion;
T = 1/Fs;             % Sampling Time
t = -2:T:2;
x = 2*tripuls(t,4);
%nfft = Fs/resulotion;
X = fft(x,N);
X = abs(X);
X = fftshift(X);
X = X/Fs;
f = -Fs/2:resulotion:Fs/2-resulotion;
figure;
plot(t,x);
title('Traingular Pulse Signal');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
actual_y_transform = 4*(sinc(2*f).^2);
%2. Derive an analytical expression for its Fourier transform. 
%3. Use Octave to calculate the Fourier transform of the signal with sampling frequency fs =10 
%Hz and resolution equal to 0.01 Hz, and then plot it together with the analytical expression on 
%one graph.
figure;
subplot(2,1,1);
plot(f,actual_y_transform);
title('Fourier Transform of the Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
legend("Analytical");
subplot(2,1,2);
plot(f,X);
title(sprintf('Fast Fourier Transform of the Signal with Fs = %d Hz',Fs));
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
legend("FFT Matlab");
%-----------------------------------------------------------
%4. Estimate the BW defined as the frequency band after which the power spectrum of the signal 
%drops to 5% of its maximum value. 
resulotion_power = 1e-5;
Fs = 10;
N_power=Fs/resulotion_power;
X_1 = fft(x,N_power);
X_1 = fftshift(X_1);
X_1 = abs(X_1);
X_1 = X_1/Fs;
f_power = linspace(-Fs/2,Fs/2,N_power);
figure;
X_squared = (X_1.^2);
X_squared_p = X_squared(f_power>=0);
f_p = f_power(f_power>=0);
plot(f_p,X_squared_p);
title('Power Spectrum of the Signal');
xlabel('Frequency (Hz)');
ylabel('X(f)^2');
xlim([0 Fs/2]);
hold on;
ind=1;
for i=1:length(X_squared_p)
    if X_squared_p(1,i) <= 0.05 * max(X_squared_p)
        ind = i;
        break
    end
end
fprintf('Max Power is  %d \n',max(X_squared_p));
fprintf('at F = %d , The power is %d \n',f_p(1,ind),X_squared_p(1,ind));
fprintf('required BW is  %d \n',f_p(1,ind));
xline(f_p(1,ind));
legend('Power Spectrum','5% of power');
hold off;
%-----------------------------------------------------------
%5. If this signal is to pass through a perfect LPF with BW= 1Hz. Plot the output of the filter in 
%the time domain along with the input signal.
figure;
passed = lowpass(x,1,Fs);
plot(t,passed);
hold on;
xlim([-2,2]);
ylim([0,2]);
plot(t,x);
xlim([-2,2]);
ylim([0,2]);
hold off;
title('Traingular Pulse Filtered with LPF at 1 Hz BW');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend('Filtered','Original');
%-----------------------------------------------------------
%6. Repeat (5) if the LPF BW is reduced to be = 0.3 Hz.
figure;
passed = lowpass(x,0.3,Fs);
plot(t,passed);
hold on;
xlim([-2.2,2.2]);
ylim([0,2]);
plot(t,x);
xlim([-2.2,2.2]);
ylim([0,2]);
hold off;
title('Traingular Pulse Filtered with LPF at 0.3 Hz BW');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend('Filtered','Original');
%-----------------------------------------------------------
%7. If this signal is to pass through a filter whose impulse response is rectangular pulse. 
%Plot the output of the filter in time domain using two methods: 
%a) Using convolution method.
figure;
t_conv = -4:1/Fs:4;
h_t = @(t) 1.*(t > -1 & t<1) + (0).*(t<=-1 & t>=1);
y_conv = conv(x,h_t(t));
y_conv = y_conv/Fs;
plot(t_conv,y_conv);
title('Filter Output using Convolution method');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend("conv Matlab");
ylim([0,3.5]);
%-----------------------------------------------------------
%b) Using multiplication in frequency domain then inverse 
%Fourier transform.
figure;
rectangular_FT = 2*sinc(2*f);
X_f_h_f = X .* rectangular_FT;
IFT = ifft(fftshift(X_f_h_f));
IFT = ifftshift(IFT);
IFT = abs(IFT)*Fs;
half_L_t_conv = (length(t_conv)-1)/2;
half_L_IFT = (length(IFT)-1)/2;
IFT = IFT(IFT >0.00001);
plot(t_conv,IFT);
title('Filter Output using Multiplication method');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend("IFFT Matlab");
%-----------------------------------------------------------
%8. Derive an analytical expression and compare it with the 
%results of (7). 
figure;
subplot(3,1,1);
plot(t_conv,y_conv);
title('Filter Output using Convolution method');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend("conv Matlab");
ylim([0,3.5]);
subplot(3,1,2);
y_conv_actual = @(t) (0.5.*t.^2+3.*t+4.5).*(t>=-3 & t<-1) + (-1.*t.^2+3).*(t>=-1 & t<1) + (0.5.*t.^2-3.*t+4.5).*(t>=1 & t<=3);  
y_convoluted = y_conv_actual(t_conv);
plot(t_conv,y_convoluted);
title('Filter Output using Convolution method');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend("Analytical");
ylim([0,3.5]);
subplot(3,1,3);
plot(t_conv,IFT);
title('Filter Output using Multiplication method');
xlabel('Time (sec)');
ylabel('Amplitude (V)');
legend("IFFT Matlab");
ylim([0,3.5]);
%-----------------------------------------------------------