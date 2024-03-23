
classdef PulseParameters
    methods (Static)
        
        % Method for calculating peak amplitude
        function peak_amplitude = calculatePeakAmplitude(psi)
            %calculate the peak amplitude
            peak_amplitude = max(abs(psi));
        end
        

        % Method for calculating temporal pulse position
        function temporal_pulse_position = calculateTemporalPulsePosition(psi, t)
            %calculate the temporal intensity
            psi_temporal_intensity = abs(psi).^2;
            %calcualte the temporal pulse position
            temporal_pulse_position = sum(t.*psi_temporal_intensity)/sum(psi_temporal_intensity); 
        end


        % Method for calculating full-width at half-maximum (FWHM)
        function FWHM = calculateFWHM(psi, t)
            %calculate the temporal intensity
            psi_temporal_intensity = abs(psi).^2;

            %full-width at half-maximum (FWHM) pulse width (not working)
            peak_intensity = max(psi_temporal_intensity); 
            half_peak_intensity = peak_intensity/2; 

            %find the position index of the left-side and right-side position with just below peak intensity value
            peak_intensity_index = find(psi_temporal_intensity==peak_intensity);
            %find the first index position (from peak) of the left-side where the intensity goes below half peak intensity
            left_index = find(psi_temporal_intensity(1:peak_intensity_index) < half_peak_intensity, 1, 'last');
            %find the first index position of the right-side where the intensity goes below half peak intensity
            right_index = peak_intensity_index + find(psi_temporal_intensity(peak_intensity_index:end) < half_peak_intensity, 1, 'first') - 1;
           
            %compute the left and right side average (using ps) 
            left_time_avg = (t(left_index) + t(left_index+1))/2;
            right_time_avg = (t(right_index-1) + t(right_index))/2;
            
            %find the FWHM
            FWHM = abs(left_time_avg - right_time_avg);
        end


        % Method for calculating frequency shift
        function frequency_shift = calculateFrequencyShift(psi, omega)
            %calculate frequency shift
            frequency_shift = trapz(omega .* (abs(fftshift(fft(psi))).^2)) / trapz((abs(fftshift(fft(psi))).^2));
        end
 
        
        % Method for calculating phase angle
        function phase_angle = calculatePhaseAngle(psi)
            %calculate phase angle
            phase_angle = mean(unwrap(angle(psi)));
        end


        % Method for calculating the chirp
        function chirp = calculateChirp(psi, t, omega, dt)
            %calculate the temporal intensity
            psi_temporal_intensity = abs(psi).^2;
            %calculate the denominator value 
            Cdenominator = sum((t.^2) .* psi_temporal_intensity) * dt;

            %calculate the numerator value 
            psi_t = ifft(fftshift(1i * omega .* fftshift(fft(psi))));
            Cnumerator = -imag(sum(t .* psi .* conj(psi_t))) * dt;

            %calculate chrip 
            chirp = Cnumerator / Cdenominator;
        end

        % Method to format the value for displaying 
        function formatted_val = formatValue(val)
            if val < 10^(-12) 
                formatted_val = '0';
            else
                if val < 10^(-3)
                    %create formatted string with scientifiy notation
                    formatted_val = sprintf('%0.4e', val);
                else
                    %create formatted string
                    formatted_val = sprintf('%0.4f', val);
                end
            end
        end
    end
end


