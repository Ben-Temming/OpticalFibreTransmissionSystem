clearvars; close('all'); clc;

L     = 28845; % Length of fiber [m]
gamma = 0.00199; % Nonlinear parameter [1/m/W]
beta2 = -2.1678e-26; % Groupe velocity dispersion parameter [sec^2/m]
loss  = 0; %0.0002; % Attenuation constant [dB/m]
alpha = loss/4.343;           % Attenuation constant [1/m]


T0   = 2.5e-11;  % Pulse width [sec]
Tm   =   25*T0;               % Time window [sec]
N    =   2^9;                % Number of modes in Fourier space
dt   =   2*Tm/N;              % Time resolution [sec]
dw   =   pi/Tm;               % Frequency resolution [rad/sec]
T    = -Tm:dt:Tm-dt;          % Time range [sec]
w    = -(pi/dt):dw:(pi/dt)-dw;% Frequency range [rad/sec]

h    = 28.845;                   % Space resolution [m]
M    = round(L/h);            % Number of space points

disp(['Updated number of samples: ' num2str(N)]);
disp(['Updated dt (time sampling rate): ' num2str(dt)]);


P0   = 0.017424;  % Peak power [Watt]
C    = 0;  % Chirp parameter
shape = "Hyperbolic secant pulse"; % Pulse shape

switch shape
    case "Gaussian pulse"
        A0 = sqrt(P0)*exp(-0.5*(1+1i*C)*(T/T0).^2);
    case "Hyperbolic secant pulse"
        A0 = sqrt(P0)*sech(T/T0).*exp(0.5i*C*(T/T0).^2);
    case "Super-Gaussian pulse"
        m = 1;  % Degree of edge sharpness
        A0 = sqrt(P0)*exp(-0.5*(1+1i*C)*(T/T0).^(2*m));
end


A       = zeros(M, N);  % Field A(z,T) : Matrix with all the calculated results
A(1,:)  = A0;           % Initial value A(0,T)

D = -alpha/2 + 0.5i*beta2*fftshift(w).^2;
%disp(D); 

for m = 2:M
    u = A(m-1,:);
    N = 1i*gamma*abs(u).^2;
    temp = fft( exp(h/2*D).*ifft(u)  );
    temp = exp(h*N).*temp;
    A(m,:) = fft(  exp(h/2*D).*ifft(temp)  );
end


psi_temporal_intensity = abs(A(end, :)).^2;  

% Plot the signal
figure;
plot(T, psi_temporal_intensity, 'LineWidth', 2);
xlabel('Time [ps]');
%ylabel('\psi(t)');
ylabel('\psi(z = 0, t)|^2 [W]');
title('Plot of \psi(t)');
grid on;

