
classdef SimulationCalculator
    properties
        T, duty_cycle, lamda0, L, alphadB, D, n2, Aeff, pulse_type
    end 

    methods 
        function obj = SimulationCalculator(T, duty_cycle, lamda0, L, alphadB, D, n2, Aeff)
            %constructor to initilaize the properties 
            obj.T = T;
            obj.duty_cycle = duty_cycle; 
            obj.lamda0 = lamda0;
            obj.L = L;
            obj.alphadB = alphadB;
            obj.D = D;
            obj.n2 = n2;
            obj.Aeff = Aeff;
        end

        %function to get simulation results 
        function [psi_evoluation, t, omega, dt, z, Nz] = performSimulation(obj, pulse_type)
            %get the speed of light 
            c = physconst('LightSpeed');
            
            %% calculate initla parameters 
            f = 1/obj.T; %frequency (Hz)
            T0 = (obj.duty_cycle/100)*obj.T; %pulse width (s)
            omega0 = (2*pi*c)/obj.lamda0; %angular frequency (rad s^-1)
            alpha = obj.alphadB/4.343; % loss in (m^-1)
            beta2 = -((obj.D*(obj.lamda0^2))/(2*pi*c)); %GVD parameter (s^2/m)
            gamma = (obj.n2*omega0)/(c*obj.Aeff); %SPM parameter (W^-1 m^-1)

            %% Define the intial pulse parameters 
            psi0 = sqrt(abs(beta2)/(gamma* T0^2)); %peak amplitude calculated from the one-solition conditoin 
            
            %calculating the linear and the non-linear length 
            LD = (T0^2)/abs(beta2); %linear length (m)
            LNL = 1/(gamma*(psi0^2)); %nonlinear length (m)

            
            %% Discretisation 
            Tmax = 40*T0; %time window width 
            fmax = 40/(2*pi*T0); %frequency window width
            
            %combining the window width and the Nyqvist criteria to avoid aliasing effects
            dt = 1/(2*fmax); %time sampling rate
            df = 1/Tmax; %frequency sampling rate
            
            %preparation for Fast Fourier Transform (FFT)
            N0 = round(Tmax/dt); % number of samples
            %adjust N0 to be in the order 2^n
            N = 2^nextpow2(N0); %set N to the 2^n closest to N0
            dt = Tmax/N; %update dt so that Tmax/dt = N
            fmax = N/(2*Tmax); %update fmax so that 2fmax/df = N
            
            %calculate the time vector (t[i] = (-N/2 + i - 1))
            t = (-N/2 : N/2 - 1)*dt; %create time vector from -N/2 to N/2 -1 and scale by dt
            %calculate the frequency vector (f[i]= (-N/2 + i - 1)
            f = (-N/2 : N/2 - 1)*df; 
            %calculate the angular vector (omega[i] = (-N/2 + i -1))
            domega = 2*pi*df; 
            omega = (-N/2 : N/2 - 1)*domega; 

            %discretisation of the space z (axial propagation direction)
            dz = (1/1000)*min(LD, LNL); %make z sampling rate a fraction of the smallest LD, LNL value 
            %disp(['dz (spatial sampling rate): ' num2str(dz)]);
            Nz = round(obj.L/dz); %calculate the number of samples 
            dz = obj.L/Nz; %update dz fit the rounded range
            %create distance vector 
            z = (0:Nz)*dz;

            %% Define the initial pulse (NEED FIXING)
            if pulse_type == "Gaussian"
                %create gaussian pulse
                psi = psi0*exp(-(t.^2)/(T0^2)); 
            else
                disp("soliton")
                %create soliton pulse
                psi = psi0*sech(t/T0); 
            end 


            %%{
            %% Apply the Split-Step Fourier Method
            %initialize matrix to store pulse at each step 
            psi_evoluation = zeros(Nz, N); 
            psi_evoluation(1, :) = psi;
            
            %calculate the dispersion term (constant over distance)
            D_hat_jw = -(alpha/2) + (1i*beta2*omega.^2)/2; 
            dispersion = exp((dz*D_hat_jw)/2); %half dispresion for more accuracy
            
            %perform SSFM for each distance step (from 1 to Nz+1, as 0 is the intial pulse)
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
            
            %%}
        end
    end
end 




%{
%% Plotting result 
%Create 3d plot for only a subset of z points to show propagation 
selected_z_indices = round(linspace(1, Nz, 20)); %select 20 points for plotting 

%selected the relevant values 
selected_psi_vals = psi_evoluation(selected_z_indices, :); 
selected_z_values = z(selected_z_indices); 

%convert z to km instead of meters 
selected_z_values_km = selected_z_values/1000; 

%calculate intensities 
selected_intensities = abs(selected_psi_vals).^2;


%create the plot
figure('Position', [100, 100, 800, 400]); % _, _, ,width, height
axes()
hold on
for i = 1:numel(selected_z_indices)
    %plot the signal for the corresponding z value
    plot3(t_ps, selected_z_values_km(i)*ones(size(t)), selected_intensities(i, :), 'Color', 'blue', 'LineWidth', 1); 
end
grid on
%add labels
xlabel('t [ps]');
ylabel('z [km]');
zlabel('\psi(z = 0, t)|^2 [W]');

%adjust rotation
view(45, 30)

%adjust axis
ylim([0, selected_z_values_km(end)]);
xticks(min(t_ps):100:max(t_ps));

%}

