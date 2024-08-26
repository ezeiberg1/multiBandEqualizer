% * Author:             Emily Zeiberg, Gabe Herman, Alvin James Bacani
% * Class:              ESE 351
% * Date:               Created 2/24/23, Last Edited 3/10/23

%% Band 1
%close all;
frequencies6 = 2*logspace(1,4,7);
frequencies6(5) = 2500;
tau1 = 1/(2*pi*frequencies6(2));
dt = 1/44100;
f = 100;
w = 2*pi*f;
t = 0:dt:8*tau1;
input = exp(1i*w*t);
x = [1 zeros(1,round(8*tau1/dt)-1)];
f_range = logspace(1,4,12);
w_range = 2*pi.*f_range;
midpoints = [41.5, 131.5 416 1566 4412 13162];
v_in = sin(2*pi*10000*t);

a_C1_ct = [1, 1/tau1];
b_C1_ct = 1/tau1;
sys1 = tf(b_C1_ct,a_C1_ct);
y_band1 = lsim(sys1, x, t);


figure(1)
subplot(2,3,1);
plot(t, y_band1);
title("Band 1 Impulse Response");
hold on;
plot(t, x);
legend(["IR", "Imp."]);
xlabel("Time [s]");
ylabel("Voltage [V]");

a_C1 = [1, -1+dt/tau1];
b_C1 = [0, dt/tau1];

H1 = zeros(1, 12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y_loop = filter(b_C1, a_C1,input_loop);
    steady_state_loop = y_loop(length(y_loop));
    H1(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(H1)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(H1)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

y_filter = filter(b_C1, a_C1, input);

figure();
plot(t, abs(y_filter));
title("Output of Band 1");
xlabel("Time [s]")
ylabel("Voltage [V]");
%% Band 2
%close all;
tau2 = 1/(2*pi*frequencies6(3));
a_C2_ct = [1, 1/tau2];
b_C2_ct = 1/tau2;
sysC2_ct = tf(b_C2_ct,a_C2_ct);

a_C2 = [1, -1+dt/tau2];
b_C2 = [0, dt/tau2];
sys_C2 = tf(b_C2,a_C2);

a_R2_ct = [1 1/tau1];
b_R2_ct = [1 0];
sysR2_ct = tf(b_R2_ct,a_R2_ct);

a_R2 = [1, dt/tau2-1];
b_R2 = [1, -1];
sys_R2 = tf(b_R2,a_R2);


y_band2 = lsim(sysR2_ct, lsim(sysC2_ct, x, t), t);

figure(1)
subplot(2,3,2);
plot(t, y_band2);
hold on;
plot(t, x);
title("Band 2 Impulse Response");
legend(["IR", "Imp"]);
xlabel("Time [s]");
ylabel("Voltage [V]");


H2 = zeros(1, 12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y_loop = filter(b_R2, a_R2, filter(b_C2, a_C2, input_loop));
    steady_state_loop = y_loop(length(y_loop));
    H2(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(H2)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Band 2");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(H2)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

y_filter = filter(b_C2, a_C2, input);
y_filter = filter(b_R2, a_R2, y_filter);

figure();
plot(t, abs(y_filter));
title("Output of Band 2");
xlabel("Time [s]")
ylabel("Voltage [V]");
%% Band 3
%close all;
tau3 = 1/(2*pi*frequencies6(4));
a_C3_ct = [1, 1/tau3];
b_C3_ct = 1/tau3;
sysC3_ct = tf(b_C3_ct,a_C3_ct);

a_C3 = [1, -1+dt/tau3];
b_C3 = [0, dt/tau3];

a_R3_ct = [1 1/tau2];
b_R3_ct = [1 0];
sysR3_ct = tf(b_R3_ct,a_R3_ct);

a_R3 = [1, dt/tau3-1];
b_R3 = [1, -1];

y_band3 = lsim(sysR3_ct, lsim(sysC3_ct, x, t), t);

figure(1)
subplot(2,3,3);
plot(t, y_band3);
hold on;
plot(t, x);
title("Band 3 Impulse Response");
legend(["IR", "Imp"]);
xlabel("Time [s]");
ylabel("Voltage [V]");

H3 = zeros(1, 12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y_loop = filter(b_R3, a_R3, filter(b_C3, a_C3, input_loop));
    steady_state_loop = y_loop(length(y_loop));
    H3(n) = steady_state_loop/input_loop(length(input_loop));
  
end

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(H3)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Band 3");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(H3)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

y_filter = filter(b_C3, a_C3, input);
y_filter = filter(b_R3, a_R3, y_filter);


figure();
plot(t, abs(y_filter));
title("Output of Band 3");
xlabel("Time [s]")
ylabel("Voltage [V]");
%% Band 4
%close all;
tau4 = 1/(2*pi*frequencies6(5));
a_C4_ct = [1, 1/tau4];
b_C4_ct = 1/tau4;
sysC4_ct = tf(b_C4_ct,a_C4_ct);

a_C4 = [1, -1+dt/tau4];
b_C4 = [0, dt/tau4];

a_R4_ct = [1 1/tau3];
b_R4_ct = [1 0];
sysR4_ct = tf(b_R4_ct,a_R4_ct);

a_R4 = [1, dt/tau4-1];
b_R4 = [1, -1];


y_band4 = lsim(sysR4_ct, lsim(sysC4_ct, x, t), t);

figure(1)
subplot(2,3,4);
plot(t, y_band4);
hold on;
plot(t, x);
title("Band 4 Impulse Response");
legend(["IR", "Imp."]);
xlabel("Time [s]");
ylabel("Voltage [V]");

H4 = zeros(1, 12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y_loop = filter(b_R4, a_R4, (filter(b_C4, a_C4, input_loop)));
    steady_state_loop = y_loop(length(y_loop));
    H4(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(H4)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Band 4");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(H4)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

y_filter = filter(b_C4, a_C4, input);
y_filter = filter(b_R4, a_R4, y_filter);


figure();
plot(t, abs(y_filter));
title("Output of Band 4");
xlabel("Time [s]")
ylabel("Voltage [V]");
%% Band 5
%close all;
tau5 = 1/(2*pi*frequencies6(6));

a_R5_ct = [1 1/tau4];
b_R5_ct = [1 0];
sysR5_ct = tf(b_R5_ct,a_R5_ct);

a_R5 = [1, dt/tau5-1];
b_R5 = [1, -1];

y_band5 = lsim(sysR5_ct, x, t);

figure(1)
subplot(2,3,5);
plot(t, y_band5);
hold on;
plot(t, x);
title("Band 5 Impulse Response");
legend(["IR", "Imp."]);
xlabel("Time [s]");
ylabel("Voltage [V]");

H5 = zeros(1, 12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y_loop = filter(b_R4, a_R4, input_loop);
    steady_state_loop = y_loop(length(y_loop));
    H5(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(H5)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Band 5");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(H5)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

y_filter = filter(b_R5, a_R5, input);

figure();
plot(t, abs(y_filter));
title("Output of Band 5");
xlabel("Time [s]")
ylabel("Voltage [V]");

figure();
plot(f_range, 20*log10(abs(H1)));
hold on;
plot(f_range, 20*log10(abs(H2)));
plot(f_range, 20*log10(abs(H3)));
plot(f_range, 20*log10(abs(H4)));
plot(f_range, 20*log10(abs(H5)));
legend(["Band 1", "Band 2", "Band 3", "Band 4",  "Band 5"]);
title("Magnitude of Frequency Response for all Bands");
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
set(gca, 'XScale', 'log')

%% Bass Boost
[inputBass, fs] = audioread("Giant Steps Bass Cut.wav");
%[inputBass, fs] = audioread("Space Station - Treble Cut.wav");
gains_bass = [30, 5, 0, 0, 0];

%Impulse Response

h_Bass = gains_bass(1)*y_band1 + gains_bass(2)*y_band2 + gains_bass(3)*y_band3 + gains_bass(4)*y_band4 + gains_bass(5)*y_band5;
figure();
plot(t, h_Bass);
hold on;
plot(t, x);
title("Bass Boost Impulse Response");
legend(["Impulse Response", "Impulse"]);
xlabel("Time [s]");
ylabel("Voltage [V]");

%Filtering
band1out = filter(b_C1, a_C1,inputBass)*gains_bass(1);
band2out = filter(b_R2, a_R2, filter(b_C2, a_C2, inputBass))*gains_bass(2);
band3out = filter(b_R3, a_R3, (filter(b_C3, a_C3, inputBass)))*gains_bass(3);
band4out = filter(b_R4, a_R4, (filter(b_C4, a_C4, inputBass))) * gains_bass(4);
band5out = filter(b_R5, a_R5, inputBass)*gains_bass(5);


sound(band1out+band2out+band3out+band4out+band5out, fs);
HBass = zeros(1,12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y1 = filter(b_C1, a_C1,input_loop)*gains_bass(1);
    y2 = filter(b_R2, a_R2,(filter(b_C2, a_C2, input_loop)))*gains_bass(2);
    y3 = filter(b_R3, a_R3, (filter(b_C3, a_C3, input_loop)))*gains_bass(3);
    y4 = filter(b_R4, a_R4, (filter(b_C4, a_C4, input_loop)))*gains_bass(4);
    y5 = filter(b_R5, a_R5, input_loop)*gains_bass(5);
    yBass = y1+y2+y3+y4+y5;
    steady_state_loop = yBass(length(yBass));
    HBass(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
plot(t, abs(yBass));
title("Output of Bass Boost");
xlabel("Time [s]")
ylabel("Voltage [V]");

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(HBass)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Bass Boost");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(HBass)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

%% Treble Boost
%[inputTreble, fs] = audioread("Giant Steps Bass Cut.wav");
[inputTreble, fs] = audioread("Space Station - Treble Cut.wav");
gains_treble = [-2, -1, 1, 1, 20];

%Impulse Response
h_treble = gains_treble(1)* y_band1 + gains_treble(2)* y_band2 + gains_treble(3) * y_band3 + gains_treble(4) * y_band4 + gains_treble(5) * y_band5;
figure();
plot(t, h_treble);
hold on;
plot(t, x);
title("Treble Boost Impulse Response");
legend(["Impulse Response", "Impulse"]);
xlabel("Time [s]");
ylabel("Voltage [V]");

%Filtering
band1out = filter(b_C1, a_C1,inputTreble) * gains_treble(1);
band2out = filter(b_R2, a_R2, filter(b_C2, a_C2, inputTreble)) * gains_treble(2);
band3out = filter(b_R3, a_R3, (filter(b_C3, a_C3, inputTreble))) * gains_treble(3);
band4out = filter(b_R4, a_R4, (filter(b_C4, a_C4, inputTreble))) * gains_treble(4);
band5out = filter(b_R5, a_R5, inputTreble)* gains_treble(5);

HTreble = zeros(1,12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y1 = filter(b_C1, a_C1,input_loop) * gains_treble(1);
    y2 = filter(b_R2, a_R2, filter(b_C2, a_C2, input_loop)) * gains_treble(2);
    y3 = filter(b_R3, a_R3, (filter(b_C3, a_C3, input_loop))) * gains_treble(3);
    y4 = filter(b_R4, a_R4, (filter(b_C4, a_C4, input_loop)))* gains_treble(4);
    y5 = filter(b_R5, a_R5, input_loop) * gains_treble(5);
    yTreble = y1+y2+y3+y4+y5;
    steady_state_loop = yTreble(length(yTreble));
    HTreble(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
plot(t, abs(yTreble));
title("Output of Treble Boost");
xlabel("Time [s]")
ylabel("Voltage [V]");

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(HTreble)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Treble Boost");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(HTreble)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')


sound(band1out+band2out+band3out+band4out+band5out, fs);
%% Unity Boost
close all;
%[inputUnity, fs] = audioread("Giant Steps Bass Cut.wav");
 [inputUnity, fs] = audioread("Space Station - Treble Cut.wav");
gains_unity = [-3, -2, 3, 3, 1]; 

%Impulse Response
h_unity = gains_unity(1) * y_band1 + gains_unity(2) * y_band2 + gains_unity(3)* y_band3 + gains_unity(4) * y_band4 + gains_unity(5) * y_band5;
figure();
plot(t, h_unity);
hold on;
plot(t, x);
title("Unity Setting Impulse Response");
legend(["Impulse Response", "Impulse"]);
xlabel("Time [s]");
ylabel("Voltage [V]");



%filtering Signal
band1out = filter(b_C1, a_C1,inputUnity) * gains_unity(1);
band2out = filter(b_R2, a_R2, (filter(b_C2, a_C2, inputUnity))) * gains_unity(2);
band3out = filter(b_R3, a_R3, (filter(b_C3, a_C3, inputUnity))) * gains_unity(3);
band4out = filter(b_R4, a_R4, (filter(b_C4, a_C4, inputUnity))) * gains_unity(4);
band5out = filter(b_R5, a_R5, inputUnity) * gains_unity(5);

unityout = band1out+band2out+band3out+band4out+band5out;

HUnity = zeros(1, 12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y1 = filter(b_C1, a_C1,input_loop) * gains_unity(1);
    y2 = filter(b_R2, a_R2, filter(b_C2, a_C2, input_loop))* gains_unity(2);
    y3 = filter(b_R3, a_R3, (filter(b_C3, a_C3, input_loop))) * gains_unity(3);
    y4 = filter(b_R4, a_R4, (filter(b_C4, a_C4, input_loop)))*gains_unity(4);
    y5 = filter(b_R5, a_R5, input_loop)*gains_unity(5);
    yunity = y1+y2+y3+y4+y5;
    steady_state_loop = yunity(length(yunity));
    HUnity(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
plot(t, abs(yunity));
title("Output of Unity Setting");
xlabel("Time [s]")
ylabel("Voltage [V]");

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(HUnity)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Unity Setting");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(HUnity)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

sound(unityout, fs);

%% Filter Background Noise
close all;
[inputNoise, fs] = audioread("Blue in Green with Siren.wav");
%sound(inputNoise, fs);
gains_noise = [5, 0, 0, 0, 0];

%Impulse Response
h_noise = gains_noise(1) * y_band1;
figure();
plot(t, h_noise);
hold on;
plot(t, x);
title("Filter for Background Noise Impulse Response");
legend(["Impulse Response", "Impulse"]);
xlabel("Time [s]");
ylabel("Voltage [V]");

%Filtering
band1out = filter(b_C1, a_C1,inputNoise)* gains_noise(1);
band2out = filter(b_R2, a_R2, (filter(b_C2, a_C2, inputNoise)))*gains_noise(2);
band3out = filter(b_R3, a_R3, (filter(b_C3, a_C3, inputNoise)))*gains_noise(3);
band4out = filter(b_R4, a_R4, (filter(b_C4, a_C4, inputNoise)))*gains_noise(4);
band5out = filter(b_R5, a_R5, inputNoise)*gains_noise(5);


sound(band1out+band2out+band3out+band4out+band5out, fs);


HNoise = zeros(1,12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y1 = filter(b_C1, a_C1,input_loop)*gains_noise(1);
    y2 = filter(b_R2, a_R2, filter(b_C2, a_C2, input_loop))*gains_noise(2);
    y3 = filter(b_R3, a_R3, (filter(b_C3, a_C3, input_loop)))*gains_noise(3);
    y4 = filter(b_R4, a_R4, (filter(b_C4, a_C4, input_loop)))*gains_noise(4);
    y5 = filter(b_R5, a_R5, input_loop)*gains_noise(5);
    yNoise = y1+y2+y3+y4+y5;
    steady_state_loop = yNoise(length(yNoise));
    HNoise(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
plot(t, abs(yNoise));
title("Output of Background Noise Filter");
xlabel("Time [s]")
ylabel("Voltage [V]");

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(HNoise)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response for Background Noise Filter");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(HNoise)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')

%% Apollo 11 Audio
close all;
[inputapollo, fs] = audioread("One Small Step.wav");
%sound(inputapollo, fs);
gains_apollo = [0, 5, .5, 0, 0];

%Impulse Response
h_apollo = gains_apollo(2)*y_band2 + gains_apollo(3) * y_band3;
figure();
plot(t, h_apollo);
hold on;
plot(t, x);
title("Apollo 11 Filter Impulse Response");
legend(["Impulse Response", "Impulse"]);
xlabel("Time [s]");
ylabel("Voltage [V]");

%Filtering
band1out = filter(b_C1, a_C1,inputapollo)*gains_apollo(1);
band2out = filter(b_R2, a_R2, filter(b_C2, a_C2, inputapollo))*gains_apollo(2);
band3out = filter(b_R3, a_R3, (filter(b_C3, a_C3, inputapollo)))*gains_apollo(3);
band4out = filter(b_R4, a_R4, (filter(b_C4, a_C4, inputapollo)))*gains_apollo(4);
band5out = filter(b_R5, a_R5, inputapollo)*gains_apollo(5);

HApollo = zeros(1,12);
for n=1:12
    w_loop = w_range(n);
    input_loop = exp(1i*w_loop*t);
    y1 = filter(b_C1, a_C1,input_loop)*gains_apollo(1);
    y2 = filter(b_R2, a_R2, filter(b_C2, a_C2, input_loop))*gains_apollo(2);
    y3 = filter(b_R3, a_R3, (filter(b_C3, a_C3, input_loop)))*gains_apollo(3);
    y4 = filter(b_R4, a_R4, (filter(b_C4, a_C4, input_loop)))*gains_apollo(4);
    y5 = filter(b_R5, a_R5, input_loop)*gains_apollo(5);
    yApollo = y1+y2+y3+y4+y5;
    steady_state_loop = yApollo(length(yApollo));
    HApollo(n) = steady_state_loop/input_loop(length(input_loop));
    
end

figure();
plot(t, abs(yApollo));
title("Output of Apollo Audio Filter");
xlabel("Time [s]")
ylabel("Voltage [V]");

figure();
subplot(2,1,1);
plot(f_range, 20*log10(abs(HApollo)));
xlabel("Frequency [Hz]");
ylabel("20log_{10}(H(\omega)) [dB]");
title("Magnitude of Frequency Response of Apollo 11 Filter");
set(gca, 'XScale', 'log')
subplot(2,1,2);
plot(f_range, angle(HApollo)/pi);
xlabel("Frequency [Hz]");
ylabel("\angle H(\omega)/\pi");
title("Phase of Frequency Response");
set(gca, 'XScale', 'log')


apollo_out = band1out+band2out+band3out+band4out+band5out;
sound(inputapollo, fs);
%audiowrite('apolloBefore.wav',inputapollo,fs);

