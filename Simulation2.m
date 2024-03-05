clc; clear; close all;

%{
All calculations are performend in SI Units (m, s, ...), so the 
values have to be converted before the calculations and be converted 
back before displaying
%}

%% Constants 
c = 2.9986*10^8; %speed of light in vacum (m/s) 
disp(['c (m/s): ' num2str(c)]);

%% define the inputs
T = 100*10^(-12); %period in seconds (from pico seconds, ps)
duty_cycle = 25; % duty cycle in % 
lamda0 = 1.55*10^(-6); %central wavelength in meter (from 1.55micro meters) 
L = 28.845*10^3; %fibre length in meter (from km) 
alphadB = 0.2*10^(-3); %loss in dB/m (from dB/km)
D = 17*10^(-6); %dispersion coefficient in s/m^2 (from ps/nm/km) 
n2 = 2.7*10^(-20); % non linear coefficient in m^2/W (from km^2/W) 
Aeff = 55*10^(-12); % effective core area in m^2 (from micro m^2) 


%% calculated inputs 
f = 1/T; %frequency (Hz)
T0 = (duty_cycle/100)*T; %pulse width (s)
omega0 = (2*pi*c)/lamda0; %angular frequency (rad s^-1)
alpha = alphadB/4.343; % loss in (m^-1)
beta2 = -((D*(lamda0^2))/(2*pi*c)); %GVD parameter (s^2/m)
gamma = (n2*omega0)/(c*Aeff); %SPM parameter (W^-1 m^-1)

%% display variables
disp(['T (s): ' num2str(T)]);
disp(['duty_cycle (%): ' num2str(duty_cycle)]);
disp(['T0 (s): ' num2str(T0)]);
disp(['lamda0 (m): ' num2str(lamda0)]);
disp(['L (m): ' num2str(L)]);
disp(['alphadB (dB/m): ' num2str(alphadB)]);
disp(['D (s/m^2): ' num2str(D)]);
disp(['n2 (m^2/W): ' num2str(n2)]);
disp(['Aeff (m^2): ' num2str(Aeff)]);
disp(['f (Hz): ' num2str(f)]);
disp(['omega0 (rad/s): ' num2str(omega0)]);
disp(['alpha (m^-1): ' num2str(alpha)]);
disp(['beta2 (s^2/m): ' num2str(beta2)]);
disp(['gamma (W^-1/m): ' num2str(gamma)]);

%% Initial pulse 
psi0 = sqrt(abs(beta2)/(gamma* T0^2)); %peak amplitude calculated from the one-solition conditoin 
%psi = psi0*sech(t/T0); 
disp(['psi0 (peak amplitude)(W^(1/2)): ' num2str(psi0)]);

%calculating the linear and the non-linear length 
LD = (T0^2)/abs(beta2); %linear length (m)
LNL = 1/(gamma*(psi0^2)); %nonlinear length (m)
disp(['LD (m): ' num2str(LD)]);
disp(['LNL (m): ' num2str(LNL)]);

%% Discretisation 
Tmax = 40*T0; %time window width 
fmax = 40/(2*pi*T0); %frequency window width
disp(['Tmax (s): ' num2str(Tmax)]);
disp(['fmax (Hz): ' num2str(fmax)]);

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
fmax = N/(2*Tmax); %update fmax so that 2fmax/df = N


disp(['Updated number of samples: ' num2str(N)]);
disp(['Updated dt (time sampling rate): ' num2str(dt)]);
disp(['Updated fmax : ' num2str(df)]);

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

% Plot the signal
%convert time to picos seconds (ps) for the plot 
t_ps = t/(10^(-12));

figure;
plot(t_ps, psi_temporal_intensity, 'LineWidth', 2);
xlabel('Time [ps]');
%ylabel('\psi(t)');
ylabel('\psi(z = 0, t)|^2 [W]');
title('Plot of \psi(t)');
grid on;

% Set y-axis limits
%ylim([0, 0.02]);

%discretisation of the space z (axial propagation direction)
dz = (1/1000)*min(LD, LNL); %make z sampling rate a fraction of the smallest LD, LNL value 
disp(['dz (spatial sampling rate): ' num2str(dz)]);
Nz = round(L/dz); %calculate the number of samples 
dz = L/Nz; %update dz fit the rounded range
%create distance vector 
z = (0:Nz)*dz;

disp(['dz (spatial sampling rate): ' num2str(dz)]);
disp(['Nz (number of samples): ' num2str(Nz)]);


%% Split-Step Fourier Method
alpha = 0; %for testing

%initialize matrix to store pulse at each step 
psi_evoluation = zeros(Nz, N); 
psi_evoluation(1, :) = psi;

%calculate the dispersion term (constant over distance)
D_hat_jw = -(alpha/2) + (1i*beta2*omega.^2)/2; 
dispersion = exp((dz*D_hat_jw)/2); %half dispresion for more accuracy


%perform SSFM for each distance step (from 1 to Nz+1, as 0 is the intial
%pulse)
for z_step = 2:Nz+1 
    %first-half dispersion 
    Psi = fftshift(fft(psi)); %fourier transform 
    Psi = Psi.*dispersion; %calculate dispersion 
    psi = ifft(fftshift(Psi));%inverse fourier transform 

    %full nonlinearity 
    N_hat = 1i*gamma*(abs(psi).^2); %nonlinear operator 
    psi = psi.*exp(dz*N_hat); %apply nonlinear operator 
    
    %second-half dispersion 
    Psi = fftshift(fft(psi)); %fourier transform 
    Psi = Psi.*dispersion; %calculate dispersion 
    psi = ifft(fftshift(Psi));%inverse fourier transform 

    %store the pulse at each step 
    psi_evoluation(z_step, :) = psi; 

end

%{
%disp(max(abs(psi_evoluation(end, :))))

%calculate the temporl intensity (|x|^2)
%psi_temporal_intensity = abs(psi).^2;  
psi_temporal_intensity = abs(psi_evoluation(end, :)).^2;  

% Plot the signal
figure;
plot(t_ps, psi_temporal_intensity, 'LineWidth', 2);
xlabel('Time [ps]');
%ylabel('\psi(t)');
ylabel('\psi(z = 0, t)|^2 [W]');
title('Plot of \psi(t) after');
grid on;
%}

% Create a meshgrid for z and time
[z_mesh, time_mesh] = meshgrid(z, t_ps);

intensities = abs(psi_evoluation).^2; 
intensities = intensities'; %transpose

% Create a 3D surface plot
figure;
surf(time_mesh, z_mesh, intensities , 'EdgeColor', 'interp');
xlabel('Time [ps]');
ylabel('Z');
zlabel('Signal Power');
title('Evolution of Signal Power along Z and Time');
view(45, 30)


%%{
%try 3d plot 
% Create a 3D surface plot using plot3
figure;
plot3(time_mesh(:), z_mesh(:), intensities(:), 'LineWidth', 2);
xlabel('Time [ps]');
ylabel('Z');
zlabel('Signal Power');
title('Evolution of Signal Power along Z and Time');
grid on;
view(45, 30)
%}



%Create 3d plot for only a subset of z points to show propagation 
selected_z_indices = round(linspace(1, Nz, 20)); %select 20 points for plotting 

%selected the relevant values 
selected_psi_vals = psi_evoluation(selected_z_indices, :); 
selected_z_values = z(selected_z_indices); 

%convert z to km instead of meters 
selected_z_values_km = selected_z_values/1000; 

%calculate intensities 
selected_intensities = abs(selected_psi_vals).^2;


%clf()
figure('Position', [100, 100, 800, 400]); % _, _, ,width, height
axes()
hold on
for i = 1:numel(selected_z_indices)
    %plot the signal for the corresponding z value
    plot3(t_ps, selected_z_values_km(i)*ones(size(t)), selected_intensities(i, :), 'Color', 'blue', 'LineWidth', 1); 
end
grid on
xlabel('t [ps]');
ylabel('z [km]');
zlabel('\psi(z = 0, t)|^2 [W]');
%view(12.6, 27.6)
view(45, 30)

ylim([0, selected_z_values_km(end)]);


% Adjust x-axis ticks and labels
%xticks(linspace(min(t_ps), max(t_ps), 10)); % You can adjust the number of ticks as needed
%xticklabels(cellstr(num2str(linspace(min(t_ps), max(t_ps), 10)'))); % Convert numerical labels to cell array of strings
