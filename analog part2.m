       %%Time specifications:
Fs=500;
t = 0:1/Fs:3;  %
   % seconds
   %%Sine wave:
   Fc = 9;                     % hertz
   x = cos(2*pi*Fc*t);
   % Plot the signal versus time:
  figure;
 plot(t,x);
   xlabel('time (in seconds)');
   title('Signal versus Time');
 zoom xon;
 figure;
%_____part1________________________________________________________________________________ 
 t = (-250:1:250)';
% The analytical solution of Fourier transform.
impulse = t==9;
impulse_1=t==-9;
%_____part2__________________________________________________________________________
subplot(2,1,1);
plot(t,impulse,'b')
hold on;
plot(t,impulse_1,'b')
hold off;
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Analytical Impulse.');
Y = fft(x);   
f = (0:length(Y)-1)*Fs/length(Y);
P = abs(Y);
n = length(x);                         
fshift = (-n/2:n/2-1)*(Fs/n);
yshift = fftshift(abs(Y));
subplot(2,1,2);plot(fshift,yshift/Fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Using FFT');
%_________part3__________________________________________________________________
 yi = max(yshift)*(1/100) ;
 %Getting the Band width of the signal.
[val,idx]=min(abs(yshift-yi));
 BW=abs(abs(fshift(idx))-9)*2;
 hold on
 plot(fshift(idx),yshift(idx)/Fs,'r*')
 disp(BW); 
 %__________part4___________________________________________________________________
%%Time specifications:
   
t = 0:1/Fs:3;  % Time vector 
 % The modulated signal.               
modul_sig=cos(2*pi*50*t).*(cos(2*pi*Fc*t));
%getting the Fourier transform.
Y_modu = fft(modul_sig);   
f_modu = (0:length(Y_modu)-1)*Fs/length(Y_modu);
P_modu = abs(Y_modu);
n = length(modul_sig);
figure;
fshift_modu = (-n/2:n/2-1)*(Fs/n);
yshift_modu = fftshift(abs(Y_modu));
plot(fshift_modu,abs(yshift_modu)/Fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('M(F)* carrier');
%getting the bandwidth.
    yi_modu = max(yshift_modu)*(1/100) ;
[val_modu,idx_modu]=min(abs(yshift_modu-yi_modu));
  BW_modu=abs(abs(fshift_modu(abs(idx_modu)))-50)*2;
  hold on
 plot(fshift_modu(idx_modu),yshift_modu(idx_modu)/Fs,'r*')
  disp(BW_modu);
 %_______part4 a___________________________________________________________
demodu_sig=modul_sig.*cos(2*pi*50*t);
Y_demodu = fft(demodu_sig);   
f_demodu = (0:length(Y_demodu)-1)*Fs/length(Y_demodu);
P_demodu = abs(Y_demodu);
n = length(demodu_sig);
fshift_demodu = (-n/2:n/2-1)*(Fs/n);
yshift_demodu = fftshift(abs(Y_demodu));
figure;
plot(fshift_demodu,yshift_demodu/Fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('demodulated in Frequency domain with No low pass filter');
%using Low pass Filter.
 yi_3 = -20 ;
 [val_3,idx_3]=min(abs(fshift_demodu-yi_3));
 DSB=yshift_demodu;
 DSB(1:idx_3)=0;
  yi_3 = 20 ;
 [val_3,idx_3]=min(abs(fshift_demodu-yi_3));
 DSB(idx_3:length(fshift_demodu))= 0;
 figure;
plot(fshift_demodu,DSB/Fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('demodulated in Frequency domain with low pass filter.');
% using inverse fourier transform.
filterd=ifftshift(DSB);
filterd=ifft(filterd);
figure
    plot(t,x,t,filterd)
    grid
   legend('Original','Filtered')
   xlabel('time (s)')
ylabel('Magnitude')
title('signal after coherent detector and the orginal signal.');
   %________________________________________DSB___________________________
% make USB only
yi_3 = 50;
yi_4 = 0;
yi_5=-50;
 [val_3,idx_3]=min(abs(fshift_modu-yi_3));
 [val_4,idx_4]=min(abs(fshift_modu-yi_4));
 [val_5,idx_5]=min(abs(fshift_modu-yi_5));
 %USB in postive side.
 ssb=yshift_modu;
 ssb(idx_4:idx_3)= 0;
 %USB in negative side.
 ssb(idx_5:idx_4)= 0;
 yi_6= max(abs(ssb))*1/100;
 [val_6,idx_6]=min(abs(ssb-yi_6));
 BW_SSB=abs(abs(fshift_modu(idx_6))-59)*2;
 disp(BW_SSB);
 %___________________________________________________C__________________
 
 figure;
 plot(fshift_modu,ssb/Fs)
 xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('SSB-USB in Frequency domain');
   hold on
 plot(fshift_modu(idx_6),ssb(idx_6)/Fs,'r*')
 hold off
 %_____________________________________________________________________
 cohernet_cos=cos(2*pi*50*t);
 filterd_1=ifftshift(ssb);
filterd_1=ifft(filterd_1);
cohernet_cos=cohernet_cos.*filterd_1;
Y_cohssb= fft(cohernet_cos);   
f_cohssb = (0:length(Y_cohssb)-1)*Fs/length(Y_cohssb);
P_cohssb = abs(Y_cohssb);
n = length(cohernet_cos);
fshift_cohssb = (-n/2:n/2-1)*(Fs/n);
yshift_cohssb = fftshift(P_cohssb);
% apply low pass filter.
 yi_3 = -20 ;
 [val_3,idx_3]=min(abs(fshift_cohssb-yi_3));
 SSB=yshift_cohssb;
 SSB(1:idx_3)=0;
  yi_3 = 20 ;
 [val_3,idx_3]=min(abs(fshift_cohssb-yi_3));
 SSB(idx_3:length(fshift_cohssb))= 0;
 filterd=ifftshift(SSB);
filterd=ifft(filterd);
figure
    plot(t,x,t,filterd)
    grid
   legend('Original','Filtered')
   xlabel('time (s)')
ylabel('Magnitude')
title('signal after coherent detector and the orginal signal.');
%______________________________________________________________SSB_________   
  ka=0.4;
   Fs=500;
t = 0:1/Fs:3;
AM_modu=(1+ka.*x).*cos(2*pi*50*t);
Y_AM= fft(AM_modu);   
f_AM = (0:length(Y_AM)-1)*Fs/length(Y_AM);
P_AM = abs(Y_AM);
n = length(AM_modu);
fshift_AM = (-n/2:n/2-1)*(Fs/n);
yshift_AM = fftshift(abs(Y_AM));
figure;
plot(fshift_AM,yshift_AM/Fs)
 xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Conventional AM in Frequency domain k_a=0.4');
yi_6= max(abs(yshift_AM))*1/100;
 [val_6,idx_6]=min(abs(abs(yshift_AM)-yi_6));
distance=0;
id=0;
for i=1:length(fshift_AM)
    if yshift_AM(i) >= yi_6
       if distance < abs(abs(fshift_AM(abs(i)))-50)
           distance=abs(abs(fshift_AM(abs(i)))-50);
           id=i;
       end
    end
end
%disp(distance); 
%disp(id);
  hold on
 plot(fshift_AM(id),yshift_AM(id)/Fs,'r*')
 BW_AM=distance*2;
 disp(BW_AM);
figure;
plot(t(1:l),AM_modu(1:l), 'b');
hold on;
envelope = abs(hilbert(AM_modu));
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('Zoomed envelope dectector VS s(t).k_a=0.4')
legend('s(t)','Envelope')
hold off
figure;
plot(t(1:l),x(1:l), 'b');
hold on;
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('Zoomed envelope dectector VS m(t).k_a=0.4')
legend('m(t)','Envelope')
hold off
figure;
plot(t,AM_modu, 'b');
hold on;
envelope = abs(hilbert(AM_modu));
% hilbert transform to calculate the envelope of the signal
plot(t,envelope, 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS s(t).k_a=0.4')
legend('s(t)','Envelope')
hold off
figure;
plot(t,x, 'b');
hold on;
% hilbert transform to calculate the envelope of the signal
plot(t,envelope, 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS m(t).k_a=0.4')
legend('m(t)','Envelope')
hold off
Y_envelope= fft(envelope);   
f_envelope = (0:length(Y_envelope)-1)*Fs/length(Y_envelope);
P_envelope = Y_envelope;
n = length(envelope);
fshift_envelope = (-n/2:n/2-1)*(Fs/n);
yshift_envelope = fftshift(abs(Y_envelope));
figure
   plot(fshift_envelope,yshift_envelope/Fs,fshift_AM,yshift_AM/Fs)
   grid
   xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('envelope dectector in frequency domain.k_a=0.4');
   legend('Envelope','modulated')
   
   %_________________________________________________________
   %___________________________________________________
   %______________________________________________________-
   ka=1.5;
   Fs = 500;
AM_modu=(1+ka.*x).*cos(2*pi*50*t);
Y_AM= fft(AM_modu);   
f_AM = (0:length(Y_AM)-1)*Fs/length(Y_AM);
P_AM = abs(Y_AM);
n = length(AM_modu);
fshift_AM = (-n/2:n/2-1)*(Fs/n);
yshift_AM = fftshift(abs(Y_AM));
figure;
plot(fshift_AM,yshift_AM/Fs)
 xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Conventional AM in Frequency domain k_a=1.5');
yi_6= max(abs(yshift_AM))*1/100;
 [val_6,idx_6]=min(abs(abs(yshift_AM)-yi_6));
distance=0;
id=0;
for i=1:length(fshift_AM)
    if yshift_AM(i) >= yi_6
       if distance < abs(abs(fshift_AM(abs(i)))-50)
           distance=abs(abs(fshift_AM(abs(i)))-50);
           id=i;
       end
    end
end
%disp(distance); 
%disp(id);
  hold on
 plot(fshift_AM(id),yshift_AM(id)/Fs,'r*')
 BW_AM=distance*2;
 disp(BW_AM);
 figure;
plot(t(1:l),AM_modu(1:l), 'b');
hold on;
envelope = abs(hilbert(AM_modu)); % hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title(' Zoomed envelope dectector VS s(t)k_a =1.5')
legend('s(t)','Envelope')
hold off 
figure;
plot(t(1:l),x(1:l), 'b');
hold on;
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title(' zoomed envelope dectector VS m(t) k_a=1.5.')
legend('m(t)','Envelope')
hold off
plot(t,AM_modu, 'b');
hold on;
envelope = abs(hilbert(AM_modu)); % hilbert transform to calculate the envelope of the signal
plot(t,envelope, 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS s(t)k_a =1.5')
legend('s(t)','Envelope')
hold off 
figure;
plot(t,x, 'b');
hold on;
% hilbert transform to calculate the envelope of the signal
plot(t,envelope, 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS m(t) k_a=1.5.')
legend('m(t)','Envelope')
hold off
%____________________________________________Envelope detector second
%way____________
  ka=0.4;
   Fs=500;
   l=200;
i=600;
t = 0:1/Fs:3;
AM_modu=(1+ka.*x).*cos(2*pi*50*t);
Y_AM= fft(AM_modu);   
f_AM = (0:length(Y_AM)-1)*Fs/length(Y_AM);
P_AM = abs(Y_AM);
n = length(AM_modu);
fshift_AM = (-n/2:n/2-1)*(Fs/n);
yshift_AM = fftshift(abs(Y_AM));
figure;
plot(fshift_AM,yshift_AM/Fs)
 xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Conventional AM in Frequency domain k_a=0.4');
yi_6= max(abs(yshift_AM))*1/100;
 [val_6,idx_6]=min(abs(abs(yshift_AM)-yi_6));
distance=0;
id=0;
for i=1:length(fshift_AM)
    if yshift_AM(i) >= yi_6
       if distance < abs(abs(fshift_AM(abs(i)))-50)
           distance=abs(abs(fshift_AM(abs(i)))-50);
           id=i;
       end
    end
end
%disp(distance); 
%disp(id);
  hold on
 plot(fshift_AM(id),yshift_AM(id)/Fs,'r*')
 BW_AM=distance*2;
 disp(BW_AM);
figure;
plot(t(1:l),AM_modu(1:l), 'b');
hold on;
d=(AM_modu.*cos(2*pi*50*t));
envelope=conv(d,exp(-t/(1/35)));
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS s(t)..ka=0.4')
legend('s(t)','Envelope')
hold off
figure;
plot(t(1:l),x(1:l), 'b');
hold on;
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS m(t)..ka=0.4')
legend('m(t)','Envelope')
hold off
Y_envelope= fft(envelope);   
f_envelope = (0:length(Y_envelope)-1)*Fs/length(Y_envelope);
P_envelope = Y_envelope;
n = length(envelope);
fshift_envelope = (-n/2:n/2-1)*(Fs/n);
yshift_envelope = fftshift(abs(Y_envelope));
figure
   plot(fshift_envelope,yshift_envelope/Fs,fshift_AM,yshift_AM/Fs)
   grid
   xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('envelope dectector in frequency domain.ka=0.4');
   legend('Envelope','modulated')
   
   
 %_________________________________________________________________
    ka=1.5;
   Fs=500;
   l=200;
t = 0:1/Fs:3;
AM_modu=(1+ka.*x).*cos(2*pi*50*t);
Y_AM= fft(AM_modu);   
f_AM = (0:length(Y_AM)-1)*Fs/length(Y_AM);
P_AM = abs(Y_AM);
n = length(AM_modu);
fshift_AM = (-n/2:n/2-1)*(Fs/n);
yshift_AM = fftshift(abs(Y_AM));
figure;
plot(fshift_AM,yshift_AM/Fs)
 xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Conventional AM in Frequency domain k_a=1.5');
yi_6= max(abs(yshift_AM))*1/100;
 [val_6,idx_6]=min(abs(abs(yshift_AM)-yi_6));
distance=0;
id=0;
for i=1:length(fshift_AM)
    if yshift_AM(i) >= yi_6
       if distance < abs(abs(fshift_AM(abs(i)))-50)
           distance=abs(abs(fshift_AM(abs(i)))-50);
           id=i;
       end
    end
end
%disp(distance); 
%disp(id);
  hold on
 plot(fshift_AM(id),yshift_AM(id)/Fs,'r*')
 BW_AM=distance*2;
 disp(BW_AM);
figure;
plot(t(1:l),AM_modu(1:l), 'b');
hold on;
d=(AM_modu.*cos(2*pi*50*t));
envelope=conv(d,exp(-t/(1/35)));
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS s(t) ka=1.5.')
legend('s(t)','Envelope')
hold off
figure;
plot(t(1:l),x(1:l), 'b');
hold on;
% hilbert transform to calculate the envelope of the signal
plot(t(1:l),envelope(1:l), 'r');
xlabel('time (s)')
ylabel('Magnitude')
title('envelope dectector VS m(t) ka=1.5.')
legend('m(t)','Envelope')
hold off
Y_envelope= fft(envelope);   
f_envelope = (0:length(Y_envelope)-1)*Fs/length(Y_envelope);
P_envelope = Y_envelope;
n = length(envelope);
fshift_envelope = (-n/2:n/2-1)*(Fs/n);
yshift_envelope = fftshift(abs(Y_envelope));
figure
   plot(fshift_envelope,yshift_envelope/Fs,fshift_AM,yshift_AM/Fs)
   grid
   xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('envelope dectector in frequency domain ka=1.5.');
   legend('Envelope','modulated')

%____________________________Envelope detector_____________________________%____________________________Envelope detector_____________________________