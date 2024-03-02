clc; clear; close all;

%% Constants 
c = 2.9986*10^8; %*10^(-12)*10^(-3);%speed of light in vacum (km/ps
disp(['c: ' num2str(c)]);
%% define the inputs
T = 100; %period in ps,  100*10^(-12); %period in seconds (from pico seconds, ps)
duty_cycle = 25; % duty cycle in % 
lamda0 = 1.55; %central wavelength in micro meters, 1.55*10^(-6); %central wavelength (from 1.55micro meters) 
L = 20; %fibre length in km,  20*10^3; %fibre length (from km) 
alphadB = 0.2; %loss in dB/km, 0.2*10^(-3); %loss in dB/m (from dB/km)

% need converting to SI unit
D = 17; %dispersion coefficient (ps/nm/km) 
n2 = 2.7*10^(-26); % non linear coeeficient (km^2/W) 
Aeff = 55; % effective core area (micro m^2) 


%% calculated inputs 
f = 1/T; %frequency in THz,  %frequency (Hz)
omega0 = (2*pi*c)/lamda0; %angular frequency 

T0 = (duty_cycle/100)*T; %pulse width in pico seconds
alpha = alphadB/4.343; % loss in km^-1
beta2 = -((D*(lamda0^2))/(2*pi*c)); %GVD parameter (supposed to be ps^2/km, not sure how)
gamma = (n2*omega0)/(c*Aeff); %SPM parameter (supposed to be W^-1 km^-1, not sure how)


%% display variables
disp(['T (ps): ' num2str(T)]);
disp(['duty_cycle (%): ' num2str(duty_cycle)]);
disp(['T0 (ps): ' num2str(T0)]);
disp(['lamda0 (um): ' num2str(lamda0)]);
disp(['L (km): ' num2str(L)]);
disp(['alphadB (dB/km): ' num2str(alphadB)]);
disp(['D (ps/nm/km): ' num2str(D)]);
disp(['n2 (km^2/W): ' num2str(n2)]);
disp(['Aeff (um^2): ' num2str(Aeff)]);
disp(['f (THz): ' num2str(f)]);
disp(['omega0 (rad/ps): ' num2str(omega0)]);
disp(['alpha (km^-1): ' num2str(alpha)]);
disp(['beta2 (ps^2/km): ' num2str(beta2)]);
disp(['gamma (W^-1/km): ' num2str(gamma)]);

%% Initial pulse 
psi0 = sqrt(abs(beta2)/(gamma* T0^2)); %peak amplitude calculated from the one-solition conditoin 
%psi0 = 0.1320;
%psi = psi0*sech(t/T0); 
disp(['psi0 (peak amplitude): ' num2str(psi0)]);

%calculating the linear and the non-linear length 
LD = (T0^2)/abs(beta2); %linear length 
LNL = 1/(gamma*(psi0^2)); %nonlinear length


%% Discretisation 
Tmax = 40*T0; %time window width 
fmax = 40/(2*pi*T0); %frequency window width
disp(['Tmax (ps): ' num2str(Tmax)]);
disp(['fmax (THz): ' num2str(fmax)]);

%combining the window width and the Nyqvist criteria to avoid aliasing
%effects
dt = 1/(2*fmax); %time sampling rate
df = 1/Tmax; %frequency sampling rate

disp(['dt (time sampling rate): ' num2str(dt)]);
disp(['df (frequency sampling rate): ' num2str(df)]);

%preparation for Fast Fourier Transform (FFT)
N0 = round(Tmax/dt); % number of samples
disp(['Number of samples: ' num2str(N0)]);
%adjust N0 to be in the order 2^n
N = 2^nextpow2(N0); %set N to the 2^n closest to N0
dt = Tmax/N; %update dt so that Tmax/dt = N
df = N/(2*Tmax); %update fmax so that 2fmax/df = N

disp(['Updated number of samples: ' num2str(N)]);
disp(['Updated dt (time sampling rate): ' num2str(dt)]);
disp(['Updated df (frequency sampling rate): ' num2str(df)]);

%calculate the time vector (t[i] = (-N/2 + i - 1))
t = (-N/2 : N/2 - 1)*dt; %create time vector from -N/2 to N/2 -1 and scale by dt
%calculate the frequency vector (f[i]= (-N/2 + i - 1)
f = (-N/2 : N/2 - 1)*df; 
%calculate the angular vector (omega[i] = (-N/2 + i -1))
domega = 2*pi*df; 
omega = (-N/2 : N/2 - 1)*domega; 

%crate the signal
psi = psi0*sech(t/T0);  
%calculate the temporl intensity (|x|^2)
psi_temporal_intensity = abs(psi).^2;  

% Plot the signal (fucked because of units)
figure;
plot(t, psi_temporal_intensity, 'LineWidth', 2);
xlabel('Time (ps)');
%ylabel('\psi(t)');
ylabel('\psi(z = 0, t)|^2 [W]');
title('Plot of \psi(t)');
grid on;

%discretisation of the space z (axial propagation direction)
dz = (1/1000)*min(LD, LNL); %make z sampling rate a fraction of the smallest LD, LNL value 
Nz = round(L/dz); %calculate the number of samples 
dz = L/Nz; %update dz fit the rounded range
