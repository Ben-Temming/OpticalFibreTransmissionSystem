classdef FibreSimulationGUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        OpticalFibreSimulationUIFigure  matlab.ui.Figure
        DropDown_TransmissionSpeed      matlab.ui.control.DropDown
        TransmissionspeedEditField      matlab.ui.control.NumericEditField
        TransmissionspeedEditFieldLabel  matlab.ui.control.Label
        FinalLabel                      matlab.ui.control.Label
        PhaseangleLabel_val_final       matlab.ui.control.Label
        FrequencyshiftLabel_val_final   matlab.ui.control.Label
        ChirpLabel_val_final            matlab.ui.control.Label
        FWHMLabel_val_final             matlab.ui.control.Label
        PulsepositionLabel_val_final    matlab.ui.control.Label
        PeakamplitudeLabel_val_final    matlab.ui.control.Label
        InitialLabel                    matlab.ui.control.Label
        PulsetypeDropDown               matlab.ui.control.DropDown
        PulsetypeDropDownLabel          matlab.ui.control.Label
        PhaseangleLabel_val_initial     matlab.ui.control.Label
        PhaseangleLabel                 matlab.ui.control.Label
        radLabel_2                      matlab.ui.control.Label
        FrequencyshiftLabel_val_initial  matlab.ui.control.Label
        FrequencyshiftLabel             matlab.ui.control.Label
        radpsLabel_2                    matlab.ui.control.Label
        ChirpLabel_val_initial          matlab.ui.control.Label
        ChirpLabel                      matlab.ui.control.Label
        radps2Label_2                   matlab.ui.control.Label
        FWHMLabel_val_initial           matlab.ui.control.Label
        FWHMLabel                       matlab.ui.control.Label
        psLabel_4                       matlab.ui.control.Label
        PulsepositionLabel_val_initial  matlab.ui.control.Label
        PulsepositionLabel              matlab.ui.control.Label
        psLabel_3                       matlab.ui.control.Label
        PeakamplitudeLabel_val_initial  matlab.ui.control.Label
        PeakamplitudeLabel              matlab.ui.control.Label
        leftsqrtWrightLabel             matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        ExitButton                      matlab.ui.control.Button
        ReportButton                    matlab.ui.control.Button
        SimulateButton                  matlab.ui.control.Button
        DropDown_AUnit                  matlab.ui.control.DropDown
        EffectivecoreareaA_texteffEditField  matlab.ui.control.NumericEditField
        EffectivecoreareaA_texteffLabel  matlab.ui.control.Label
        DropDown_NUnit                  matlab.ui.control.DropDown
        Nonlinearcoefficientn_2EditField  matlab.ui.control.NumericEditField
        Nonlinearcoefficientn_2Label    matlab.ui.control.Label
        DropDown_DUnit                  matlab.ui.control.DropDown
        DispersioncoefficientDEditField  matlab.ui.control.NumericEditField
        DispersioncoefficientDEditFieldLabel  matlab.ui.control.Label
        DropDown_LossUnit               matlab.ui.control.DropDown
        Lossalpha_textdBEditField       matlab.ui.control.NumericEditField
        Lossalpha_textdBLabel           matlab.ui.control.Label
        DropDown_FibreLengthUnit        matlab.ui.control.DropDown
        FibrelengthLEditField           matlab.ui.control.NumericEditField
        FibrelengthLEditFieldLabel      matlab.ui.control.Label
        DropDown_WavelengthUnit         matlab.ui.control.DropDown
        Wavelengthlambda_0EditField     matlab.ui.control.NumericEditField
        Wavelengthlambda_0EditFieldLabel  matlab.ui.control.Label
        DutycycleEditField              matlab.ui.control.NumericEditField
        DutycycleEditFieldLabel         matlab.ui.control.Label
        OutputsLabel                    matlab.ui.control.Label
        InputsLabel                     matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        %store data for report generation 
        t_ps_array %time array 
        initial_intensity_array %inital intensity arrray 
        selected_z_values_km_array %selected z values in km
        selected_intensities_array %intensities for each selected z value
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SimulateButton
        function SimulateButtonPushed(app, event)
            %disable simulation button to avoid clicking while calculating
            app.SimulateButton.Enable = "off";

            %% read all inputs
            % get the transmission speed value and check if field is empty
            t_speed_val = app.TransmissionspeedEditField.Value; 
            if isempty(t_speed_val)
                % Display an error message 
                errordlg('Transmission speed field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return; 
            end 
            %get the transmission speed unit
            t_speed_unit = app.DropDown_TransmissionSpeed.Value; 

            % get duty cycle value and check if field is empty
            duty_cycle_val = app.DutycycleEditField.Value;
            if isempty(duty_cycle_val)
                % Display an error message
                errordlg('Duty cycle field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end

            %get wavelength value and check if field is empty 
            Wavelength = app.Wavelengthlambda_0EditField.Value;
            if isempty(Wavelength)
                % Display an error message
                errordlg('Wavelength field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end
            %get the wavelength unit
            Wavelength_unit = app.DropDown_WavelengthUnit.Value; 

            % get Fibrelength value and check if field is empty
            Fibrelength = app.FibrelengthLEditField.Value;
            if isempty(Fibrelength)
                % Display an error message
                errordlg('Fibrelength field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end
            %get the fibre length unit
            Fibrelength_unit = app.DropDown_FibreLengthUnit.Value; 

            % get loss value and check if field is empty
            loss = app.Lossalpha_textdBEditField.Value;
            if isempty(loss)
                % Display an error message
                errordlg('Loss field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end
            %get the loss unit
            loss_unit = app.DropDown_LossUnit.Value; 

            % get dispersion value and check fi field is empty
            dispersion = app.DispersioncoefficientDEditField.Value;
            if isempty(dispersion)
                % Display an error message
                errordlg('Dispersion coefficient field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end
            %get the dispersion unit
            dispersion_unit = app.DropDown_DUnit.Value;

            %get non linear value and check if field is empty
            non_linear = app.Nonlinearcoefficientn_2EditField.Value;
            if isempty(non_linear)
                % Display an error message
                errordlg('Nonlinear coefficient field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end
            %get the non linear unit
            non_linear_unit = app.DropDown_NUnit.Value;  

            %get core area value and check if field is empty
            area = app.EffectivecoreareaA_texteffEditField.Value;
            if isempty(area)
                % Display an error message
                errordlg('Effective core area field cannot be empty.', 'Input Error', 'modal');
                app.SimulateButton.Enable = "on";
                return;
            end
            %get the core area unit
            area_unit = app.DropDown_AUnit.Value; 
           
            %get the pulse type
            pulse_type = app.PulsetypeDropDown.Value; 

            %% Convert all units to SI units 
            %if transmission speed is not in terms of b/s, convert
            t_speed = t_speed_val; 
            switch t_speed_unit
                case "Gb/s"
                    t_speed = t_speed*10^9; 
                case "Mb/s"
                    t_speed = t_speed*10^6; 
                case "Kb/s"
                    t_speed = t_speed*10^3; 
            end
       
            %if wavelength is not in terms of m, convert 
            lamda0 = Wavelength; 
            switch Wavelength_unit
                case "mm"
                    lamda0 = lamda0*10^(-3); 
                case "μm"
                    lamda0 = lamda0*10^(-6);
                case "nm"
                    lamda0 = lamda0*10^(-9);
            end

            %if fibre length is not in terms of m, convert 
            L = Fibrelength; 
            if Fibrelength_unit == "km"
                L  = L*10^3;
            end

            %if loss is not in terms of dB/m, convert
            alphadB = loss; 
            if loss_unit == "dB/km"
                alphadB = alphadB*10^(-3); 
            end

            %if dispersion not in s/m^2
            D = dispersion; 
            if dispersion_unit == "ps/nm/km"
                D = D*10^(-6); 
            end

            %if non linear coefficient not in terms of  m^2/W, convert 
            n2 = non_linear; 
            if non_linear_unit == "km²/W"
                n2 = n2*10^(-20); 
            end
            
            %if area not in terms of m², convert
            Aeff = area; 
            if area_unit == "μm²" 
                Aeff = Aeff*10^(-12);
            elseif area_unit == "mm²"
                Aeff = Aeff*10^(-6); 
            end
                       
            %% Perform Simulation Calculations 
            %calculate the period from the transmission speed
            T = 1/t_speed; 
            
            %setup the simulation class
            Simulation = SimulationCalculator(T, duty_cycle_val, lamda0, L, alphadB, D, n2, Aeff); 
            %get the result of the simulation 
            [psi_evoluation, t, omega, dt, z, Nz] = Simulation.performSimulation(pulse_type); 
                
            %get t in ps 
            t_ps = t/(10^(-12));
            %get dt in ps
            dt_ps = dt*(10^12);
            %get omega in ps
            omega_rad_ps = omega * 10^(-12);

            %% Calculate Initial and final parameters 
            %select initial and final pulse
            psi_initial = psi_evoluation(1, :); 
            psi_final = psi_evoluation(end, :); 

            %calculate inital and final peak amplitude 
            peak_amplitude_initial = PulseParameters.calculatePeakAmplitude(psi_initial); 
            peak_amplitude_final = PulseParameters.calculatePeakAmplitude(psi_final);

            %calculate inital and final temporal pulse position 
            temporal_pulse_position_initial = PulseParameters.calculateTemporalPulsePosition(psi_initial, t_ps);
            temporal_pulse_position_final = PulseParameters.calculateTemporalPulsePosition(psi_final, t_ps);

            %calculate the inital and final FWHM 
            FWHM_initial = PulseParameters.calculateFWHM(psi_initial, t_ps);
            FWHM_final = PulseParameters.calculateFWHM(psi_final, t_ps); 

            %chirp 
            chirp_inital = PulseParameters.calculateChirp(psi_initial, t_ps, omega_rad_ps, dt_ps); 
            chirp_final = PulseParameters.calculateChirp(psi_final, t_ps, omega_rad_ps, dt_ps); 

            %frequency shift 
            freq_shift_initial = PulseParameters.calculateFrequencyShift(psi_initial, omega); 
            freq_shift_initial = freq_shift_initial* 10^(-12); % Get it in ps
            freq_shift_final = PulseParameters.calculateFrequencyShift(psi_final, omega); 
            freq_shift_final = freq_shift_final* 10^(-12); % Get it in ps

            %phase angle 
            phase_angle_intial = PulseParameters.calculatePhaseAngle(psi_initial); 
            phase_angle_final = PulseParameters.calculatePhaseAngle(psi_final); 
            

            %% Display Parameters 
            app.PeakamplitudeLabel_val_initial.Text = PulseParameters.formatValue(peak_amplitude_initial);
            app.PeakamplitudeLabel_val_final.Text = PulseParameters.formatValue(peak_amplitude_final);
            app.PulsepositionLabel_val_initial.Text = PulseParameters.formatValue(temporal_pulse_position_initial); 
            app.PulsepositionLabel_val_final.Text = PulseParameters.formatValue(temporal_pulse_position_final); 
            app.FWHMLabel_val_initial.Text = PulseParameters.formatValue(FWHM_initial); 
            app.FWHMLabel_val_final.Text = PulseParameters.formatValue(FWHM_final); 
            app.ChirpLabel_val_initial.Text = PulseParameters.formatValue(chirp_inital); 
            app.ChirpLabel_val_final.Text = PulseParameters.formatValue(chirp_final); 
            app.FrequencyshiftLabel_val_initial.Text = PulseParameters.formatValue(freq_shift_initial); 
            app.FrequencyshiftLabel_val_final.Text = PulseParameters.formatValue(freq_shift_final);
            app.PhaseangleLabel_val_initial.Text = PulseParameters.formatValue(phase_angle_intial); 
            app.PhaseangleLabel_val_final.Text = PulseParameters.formatValue(phase_angle_final); 

            %% Create and display the plot
            %Create 3d plot for only a subset of z points to show propagation
            selected_z_indices = round(linspace(1, Nz, 20)); %select 20 points for plotting
            %selected the relevant values 
            selected_psi_vals = psi_evoluation(selected_z_indices, :); 
            selected_z_values = z(selected_z_indices); 
            %convert z to km instead of meters 
            selected_z_values_km = selected_z_values/1000;
            %calculate intensities 
            selected_intensities = abs(selected_psi_vals).^2;

            %clear the UIAxes
            cla(app.UIAxes);
            %plot on the UIAxis 
            hold(app.UIAxes, "on"); 
            for i = 1:numel(selected_z_indices)
                %plot the signal for the corresponding z value
                 plot3(app.UIAxes, t_ps, selected_z_values_km(i)*ones(size(t_ps)), selected_intensities(i, :), 'Color', 'blue', 'LineWidth', 1);
            end
            %add grid to plot
            grid(app.UIAxes, 'on');

            %add labels
            xlabel(app.UIAxes, 't [ps]');
            ylabel(app.UIAxes, 'z [km]');
            zlabel(app.UIAxes, '|\psi(z = 0, t)|^2 [W]');

            %adjust rotation
            view(app.UIAxes, 45, 30)

            %adjust axis
            ylim(app.UIAxes, [0, selected_z_values_km(end)]);
            %xticks(app.UIAxes, min(t_ps):100:max(t_ps));
            %create 10 evenly spaced ticks on the x axis
            plot_t_step = round((abs(min(t_ps) - max(t_ps))/10)/10)*10;
            xticks(app.UIAxes, min(t_ps):plot_t_step:max(t_ps));

            %add title
            title(app.UIAxes, "Signal Propagation Through Optical Fibre");

            %% Store data for report generation 
            %store t 
            app.t_ps_array = t_ps; 
            %store intial intensity
            app.initial_intensity_array = abs(psi_initial).^2;
            %store selected z 
            app.selected_z_values_km_array = selected_z_values_km;
            %store selected intensities 
            app.selected_intensities_array = selected_intensities; 

            %% Enable buttons and display plot
            app.SimulateButton.Enable = "on"; %turn simulation button back on
            % After calculations, enable the ReportButton
            app.ReportButton.Enable = 'on';
            %display the plot 
            app.UIAxes.Visible = 'on';
        end

        % Button pushed function: ReportButton
        function ReportButtonPushed(app, event)
            %% Disable all buttons
            app.SimulateButton.Enable = "off";
            app.ReportButton.Enable = "off";

            %% read the input parameters from the GUI  
            t_speed_val_str = num2str(app.TransmissionspeedEditField.Value); 
            duty_cycle_val_str = num2str(app.DutycycleEditField.Value); 
            lamda0_val_str =  num2str(app.Wavelengthlambda_0EditField.Value); 
            fibre_len_val_str =  num2str(app.FibrelengthLEditField.Value);
            alphadB_val_str = num2str(app.Lossalpha_textdBEditField.Value); 
            D_val_str = num2str(app.DispersioncoefficientDEditField.Value); 
            n2_val_str = num2str(app.Nonlinearcoefficientn_2EditField.Value); 
            Aeff_val_str = num2str(app.EffectivecoreareaA_texteffEditField.Value);    
            pulse_type_str = app.PulsetypeDropDown.Value;

            %% read the units from the GUI
            t_speed_unit = app.DropDown_TransmissionSpeed.Value; 
            Wavelength_unit = app.DropDown_WavelengthUnit.Value; 
            Fibrelength_unit = app.DropDown_FibreLengthUnit.Value; 
            loss_unit = app.DropDown_LossUnit.Value; 
            dispersion_unit = app.DropDown_DUnit.Value;
            non_linear_unit = app.DropDown_NUnit.Value; 
            area_unit = app.DropDown_AUnit.Value;

            %% create the corresponding SI unit string 
            %create the transmission speed unit string 
            switch t_speed_unit
                case "Gb/s"
                    t_speed_unit_str = "\\giga\\bit\\per\\second";
                case "Mb/s"
                    t_speed_unit_str = "\\mega\\bit\\per\\second";
                case "Kb/s"
                    t_speed_unit_str = "\\kilo\\bit\\per\\second";
                case "b/s"
                    t_speed_unit_str = "\\bit\\per\\second";
            end

            %create the wavelength unit string
            switch Wavelength_unit
                case "m"
                    lamda0_unit_str = '\\meter';
                case "mm"
                    lamda0_unit_str = '\\milli\\meter';  
                case "μm"
                    lamda0_unit_str = '\\micro\\meter'; 
                case "nm"
                    lamda0_unit_str = '\\nano\\meter'; 
            end

            %create the fibre length unit string
            if Fibrelength_unit == "km"
                fibre_len_unit_str = '\\kilo\\meter'; 
            else
                fibre_len_unit_str = '\\meter';
            end

            %create the loss unit string 
            if loss_unit == "dB/km"
                alphadB_unit_str = '\\deci\\bel\\per\\kilo\\meter';
            else
                alphadB_unit_str = '\\deci\\bel\\per\\meter';
            end

            %create dispersion unit string 
            if dispersion_unit == "ps/nm/km"
                D_unit_str = '\\pico\\second\\per\\nano\\meter\\per\\kilo\\meter';
            else
                D_unit_str = '\\second\\per\\meter\\squared';
            end
            
            %create non linear unit string 
            if non_linear_unit == "km²/W"
                n2_unit_str = '\\kilo\\meter\\squared\\per\\watt';
            else
                n2_unit_str = '\\meter\\squared\\per\\watt';
            end

            %create area unit string 
            if area_unit == "μm²" 
                Aeff_unit_str = '\\micro\\meter\\squared'; 
            elseif area_unit == "mm²"
                Aeff_unit_str = '\\milli\\meter\\squared'; 
            else
                Aeff_unit_str = '\\meter\\squared'; 
            end

            %% Get the output parameters
            %read the table parameters from the GUI 
            inital_peak_a_str = app.PeakamplitudeLabel_val_initial.Text; 
            final_peak_a_str = app.PeakamplitudeLabel_val_final.Text;
            initial_pulse_pos_str = app.PulsepositionLabel_val_initial.Text; 
            final_pulse_pos_str = app.PulsepositionLabel_val_final.Text; 
            initial_FWHM_str = app.FWHMLabel_val_initial.Text; 
            final_FWHM_str = app.FWHMLabel_val_final.Text; 
            inital_chirp_str = app.ChirpLabel_val_initial.Text; 
            final_chirp_str = app.ChirpLabel_val_final.Text; 
            initial_freq_shift_str = app.FrequencyshiftLabel_val_initial.Text;
            final_freq_shift_str = app.FrequencyshiftLabel_val_final.Text; 
            initial_phase_ang_str = app.PhaseangleLabel_val_initial.Text; 
            final_phase_ang_str = app.PhaseangleLabel_val_final.Text; 
            
            %% Specify the data for the initial pulse plot
            %select values to appropriatly scale the plot
            initial_pulse_xmin_str = num2str(min(app.t_ps_array(:))); 
            initial_pulse_xmax_str  = num2str(abs(min(app.t_ps_array(:))));
            initial_pulse_ymin_str = num2str(min(app.initial_intensity_array(:)));
            initial_pulse_ymax_str = num2str(max(app.initial_intensity_array(:))*1.1);
             
            % Get the size of the array
            [t_rows, t_cols] = size(app.t_ps_array);
            % Initialize an empty cell array to store the data
            initial_pulse_data = cell(t_cols, 1);
            % loop over the time and intensity array
            for i = 1:t_cols
                %construct the string representing time and intensity
                time_intensity_str = ['(' num2str(app.t_ps_array(i)) ', ' num2str(app.initial_intensity_array(i)) ')'];
                %store the string in the cell array
                initial_pulse_data{i} = time_intensity_str;
            end
            
            %% Generate the pulse propagation plot
            %define the empty cell array
            pulse_propagation_plot = {};
            pulse_propagation_plot = [pulse_propagation_plot; {'%%plot the pulse at each selected distance'}]; 

            %Get the size of the array
            [z_rows, z_cols] = size(app.selected_z_values_km_array);

            %loop over each selected z value and generate a plot 
            for i = 1:z_cols
                %get the z value
                z_val = app.selected_z_values_km_array(i); 

                % Initialize an empty cell array to store the plot
                plot_data = cell(t_cols+2, 1);
                %add the start of the plot 
                plot_data{1} = '\\addplot3 [area plot] coordinates {'; 
                %loop over the time and intensity array
                for j = 1:t_cols
                    %construct the string representing time, distance and intensity
                    data_str = ['    (' num2str(app.t_ps_array(j)) ', ' num2str(z_val) ', ' num2str(app.selected_intensities_array(i, j)) ')'];
                    %store the string in the cell array
                    plot_data{j+1} = data_str;
                end
                %add the end of the plot
                plot_data{end} = '};';
 
               %add the plot data to the plot 
               pulse_propagation_plot = [pulse_propagation_plot; plot_data]; 
            end

            %% Generate cell array for latex report 
            report = [{
                '%% A4 paper size, flush left equations, article type document', 
                '\\documentclass[a4paper,fleqn]{article}', 
                '%% this package is included to adjust page settings',
                '\\usepackage{geometry}', 
                '\\geometry{margin=0.85in}',
                '\\usepackage[x11names]{xcolor}',
                '\\usepackage{amsmath}',
                '\\usepackage{physics} % for partial derivatives',
                '\\usepackage{hyperref}',
                '\\hypersetup{colorlinks,',
                '    citecolor={Blue1},',
                '    filecolor={Blue1},',
                '    linkcolor={Blue1},',
                '    urlcolor={Blue1}}',
                '\\urlstyle{same}',
                '\\usepackage{fancyhdr}',
                '\\usepackage{xspace}',
                '\\usepackage{graphicx}',
                '\\usepackage{minitoc}',
                '\\usepackage{siunitx}',
                '\\usepackage{tikz}',
                '\\usepackage[RPvoltages]{circuitikz}',
                '\\usepackage{varioref}',
                '\\usepackage{wrapfig}',
                '\\usetikzlibrary{calc}',
                '\\usetikzlibrary{arrows}',
                '\\usepackage{pgfplots}',
                ''
                '\\parskip=0.15cm',
                '\\parindent=0.0cm',
                '',
                '\\pdfinfo {',
                '  /Title (Optical Fibre Transmission System Simulation)',
                '  /Author (Ben Temming)',
                '  /Subject (EE4545 Assignment 1)',
                '}',
                '',
                '\\newcommand{\\shorttitle}{}',
                '\\pagestyle{fancy}',
                '\\renewcommand{\\headrulewidth}{0.5pt}',
                '\\renewcommand{\\footrulewidth}{0.5pt}',
                '\\rhead{\\shorttitle}',
                '\\chead{EE4545}',
                
                '\\cfoot{Page \\thepage\\xspace of \\pageref{LastPage}}',
                ''
                '\\fancypagestyle{plain}{',
                '  \\chead{EE4545}',
                '  \\rhead{\\shorttitle}',
                '  \\cfoot{Page \\thepage\\xspace of \\pageref{LastPage}}',
                '}',
                '',
                '\\title{Optical Fibre Transmission System Simulation}',
                '\\author{Ben Temming}',
                '\\date{EE4545 Assignment 1}',
                '',
                '\\begin{document}',
                '',
                '\\maketitle',
                '',
                '\\begingroup',
                '\\hypersetup{linkcolor=Blue1}',
                '\\tableofcontents',
                '\\endgroup',

                '\\clearpage  %% go to next page',
                '',
                '\\section{Introduction}',
               
                ['This report details the design and simulation process of an ' ...
                'optical fibre transmission system using the Split-Step Fourier Method ' ...
                '(SSFM) to solve the Nonlinear Schrödinger Equation (NLSE). ' ...
                'The system model equation, signal and system parameters and' ...
                ' pulse parameter computations are discussed.  '],
                '', 
                '\\section{System Model Equation}',
                ['The optical pulse propagation in the fibre system is modelled' ...
                ' by the nonlinear partial differential equation \\cite{non_linear_optics}: '],
                '\\vspace{-0.045in}',
                '\\begin{align}',
                '\\label{nlse}',
                ['\\frac{\\partial\\psi}{\\partial z}+\\frac{\\alpha}{2}\\psi+\\frac' ...
                '{\\jmath\\beta_2}{2}\\frac{\\partial^2\\psi}{\\partial t^2}-\\jmath\\gamma ' ...
                '|\\psi |^2\\psi=0'],
                '\\end{align}',
                '',
                ['Using the Split-Step Fourier Method (SSFM), this equation can be solved. ' ...
                'To get to a solution, it can be assumed that the dispersion and nonlinearity ' ...
                'act independently over short distances, even though they act together ' ...
                'along the length of the fibre. Using this knowledge it can be shown that ' ...
                'the optical pulse after propagating a single distance step can be described' ...
                ' by the following equation: '],
                '\\vspace{-0.045in}',
                '\\begin{align}',
                '\\label{psi_update}',
                ['\\psi(z+dz,t) \\approx \\exp\\left(\\frac{dz}{2}\\hat{D}\\right)' ...
                '\\exp(dz\\hat{N})\\exp\\left(\\frac{dz}{2}\\hat{D}\\right)\\psi(z,t)'],
                '\\end{align}',
                '',
                'where ',
                ['$z$ is the distance travelled, $dz$ is the distance step, $\\hat{D}$ ' ...
                'is the linear operator given by'],
                '\\begin{align}',
                '\\label{linear}',
                '    \\hat{D} \\equiv -\\frac{\\alpha}{2}-\\frac{\\jmath\\beta_2}{2}\\pdv[2]{t}',
                '\\end{align}',
                'and $\\hat{N}$ is the non-linear operator given by',
                '\\begin{align}',
                '\\label{nonlinear}',
                '    \\hat{N} \\equiv \\jmath\\gamma |\\psi |^2',
                '\\end{align}',
                '',
                ['Equation \\eqref{psi_update} represents the update rule for the optical ' ...
                'pulse, incorporating both dispersion and nonlinearity effects. Implementation ' ...
                'of Equation \\eqref{psi_update} in MATLAB, followed by its iterative application,' ...
                ' allows for the modelling of the pulse evolution along the length of the fibre.'],
        
                 %write the signal and system parameters
                '\\section{Signal and System Parameters}',
                '',
                'The simulation is performed with the following input parameters: ',
                '\\begin{itemize}',                
                '\\item \\textbf{Transmission speed} = \\( \\SI{',
                t_speed_val_str,
                '    }{',
                t_speed_unit_str,
                '    } \\)',
                '    \\item \\textbf{Duty cycle} = \\( ',
                duty_cycle_val_str,
                '    \\) \\%%',
                '    \\item \\textbf{Wavelength (\\( \\lambda_0 \\))} = \\( \\SI{',

                lamda0_val_str,
                '    }{',
                lamda0_unit_str,            
                '    } \\)',
                '    \\item \\textbf{Fibre length (\\( L \\))} = \\( \\SI{',
                fibre_len_val_str,
                '    }{',
                fibre_len_unit_str,
                '    } \\)',
                '    \\item \\textbf{Loss (\\( \\alpha_{\\text{dB}} \\))} = \\( \\SI{',
                alphadB_val_str,
                '    }{',
                alphadB_unit_str,
                '    } \\)',
                '    \\item \\textbf{Dispersion coefficient (\\( D \\))} = \\( \\SI{',
                D_val_str,
                '    }{',
                D_unit_str,
                '    } \\)',
                '    \\item \\textbf{Nonlinear coefficient (\\( n_2 \\))} = \\( \\SI{',
                n2_val_str, 
                '    }{',
                n2_unit_str,
                '    } \\)',
                '    \\item \\textbf{Effective core area (\\( A_{\\text{eff}} \\))} = \\( \\SI{',
                Aeff_val_str,
                '    }{',
                Aeff_unit_str,
                '    } \\)',
                '    \\item \\textbf{Pulse type} = ', 
                pulse_type_str,
                '\\end{itemize}',
                '',
                '\\newpage',

                '\\section{Simulation Result}',
                '',
                ['This section shows the result of the simulation based on ' ...
                'the parameters defined earlier. Utilizing the parameters, ' ...
                'the initial pulse is generated, serving as the starting point of' ...
                ' the simulation. '],
                '',
                '%% initial pulse figure ',
                '\\begin{figure}[!htpb]',
                '\\centering',
                '\\begin{tikzpicture}[scale=1]',
                '    \\begin{axis}[',
                '        thick,',
                '        scale only axis,',
                '        height=5cm,',
                '        width=12cm,',
                '        grid=none,',
                '        xmin= ',
                initial_pulse_xmin_str, 
                ',',
                '        xmax=',
                initial_pulse_xmax_str, 
                ',', 
                '        ymin=',
                initial_pulse_ymin_str, 
                ',', 
                '        ymax=', 
                initial_pulse_ymax_str, 
                ',', 
                '        tick align=outside,',
                '        yticklabel style={',
                '            /pgf/number format/fixed,',
                '            /pgf/number format/precision=3,',
                '            /pgf/number format/fixed zerofill',
                '        },',
                '        scaled y ticks=false,',
                '        xlabel={$t$ [\\si{\\pico\\second}]},',
                '        ylabel={$|\\psi(z=0,t)|^2$ [\\si{\\watt}]},',
                '        xlabel style={rotate=0,yshift=0cm,xshift=0.0cm, anchor=near ticklabel},',
                '        ylabel style={rotate=0,yshift=0cm,xshift=0cm, anchor=near ticklabel},    ',
                '    ]',
                '    \\addplot [',
                '        no markers,',
                '        very thick,',
                '        blue,',
                '        smooth,',
                '    ] coordinates {',
                }; 
                %insert inital plot data here
                initial_pulse_data;
                {
                '      };',
                '    \\end{axis}',
                '\\end{tikzpicture}',
                '\\vspace{-0.1in}',
                '\\caption{Initial pulse.}',
                '\\label{initial_pulse}',
                '\\end{figure}',
                '%%description of the figure',
                ['The initial pulse, as shown in Figure \\ref{initial_pulse}, ' ...
                'represents the optical signal at the beginning of the simulation. ']
                ['This pulse is the pulse generated at the transmitter before it is ' ...
                'send though the optical fibre. '],
                '',
                '%%3d plot of pulse throughout the fibre',
                '\\vspace{0.1in}',
                '\\begin{figure}[!htpb]',
                '\\centering',
                '\\begin{tikzpicture}',
                '\\begin{axis}[', 
                '    scale only axis,', 
                '    height=8cm,', 
                '    width=14cm,', 
                '    zmin=0,', 
                '    xlabel={$t$ [\\si{\\pico\\second}]},', 
                '    ylabel={$z$ [\\si{\\kilo\\meter}]},', 
                '    zlabel={$|\\psi(z,t)|^2$ [\\si{\\watt}]},', 
                '    view = {45}{45},', 
                '    area plot/.style={', 
                '        fill opacity=0.25,', 
                '        draw=blue,thick,', 
                '        fill=blue!30!white,', 
                '        mark=none,', 
                '    }', 
                ']',
                 };
                 %insert the 3d plot here
                 pulse_propagation_plot;
                 {
                 '\\end{axis}'
                '\\end{tikzpicture}'
                '\\caption{Pulse propagation through the fibre.}'
                '\\label{pulse_prop}'
                '\\end{figure}'

                ['Figure \\ref{pulse_prop} shows the evolution of the initial ' ...
                'pulse, as shown in Figure \\ref{initial_pulse}, as it traverses ' ...
                'the optical fibre. This provides insights into the influence of ' ...
                'key parameters such as fibre length, dispersion coefficient, ' ...
                'nonlinear coefficient and loss on the pulse''s characteristics.'],
                '',
                '\\newpage',
                '',
                ['To get a better understanding of how the pulse has changed as it ' ...
                'has propagated through the fibre a few pulse parameters such as ' ...
                'the peak amplitude and the temporal pulse position are calculated ' ...
                'for the initial and the final pulse.'],
                '\\begin{table}[!htb]',
                '\\centering',
                ['\\caption{Comparison of parameters between the initial ' ...
                'pulse and the final pulse.}'],
                '\\vspace{0.05in}',
                '\\begin{tabular}{|l|c|c|}',
                '    \\hline',
                '    Parameter & Initial pulse & Final pulse \\\\ \\hline\\hline',
                '    Peak amplitude [\\si{\\sqrt\\watt}] & \\num{',
                inital_peak_a_str,
                '    } & \\num{',
                final_peak_a_str,
                '    } \\\\',
                '    Pulse position [\\si{\\pico\\second}] & $\\num{',
                initial_pulse_pos_str,
                '    }$ & $\\num{',
                final_pulse_pos_str,
                '    }$ \\\\',
                '    FWHM [\\si{\\pico\\second}] & \\num{',
                initial_FWHM_str,
                '    } & \\num{',
                final_FWHM_str,
                '    } \\\\',
                '    Chirp [\\si{\\radian/\\pico\\second\\squared}] & $\\num{',
                inital_chirp_str,
                '    }$ & \\num{',
                final_chirp_str,
                '    } \\\\',
                '    Frequency shift [\\si{\\radian/\\pico\\second}] & $\\num{',
                initial_freq_shift_str,
                '    }$ & $\\num{',
                final_freq_shift_str,
                '    }$ \\\\',
                '    Phase angle [\\si{\\radian}] & $\\num{',
                initial_phase_ang_str,
                '    }$ & $\\num{',
                final_phase_ang_str,
                '    }$ \\\\',
                '    \\hline',
                '\\end{tabular}',
                '\\label{table:parameter_comparison}',
                '\\vspace{0.05in}',
                '\\end{table}',
                '',
                ['Table \\ref{table:parameter_comparison} presents a comparison' ...
                ' of the parameters between the initial and final pulses. ' ...
                'It is important to note that values smaller than $\\num{e-12}$ ' ...
                'are approximated to $\\num{0}$ as the calculations are affected ' ...
                'by the floating point error for smaller numbers. The parameters ' ...
                'characterize the temporal and spectral properties of the pulse ' ...
                'at the beginning and the end of the optical fibre. Comparing ' ...
                'these parameters gives valuable insight into the behaviour of ' ...
                'the pulse and the effects of dispersion, attenuation and ' ...
                'nonlinearities on its characteristics. '],      
                '',
                '',
                '\\section{Conclusion}',
                ['In conclusion, this report examines the propagation of a pulse through ' ...
                'an optical fibre, simulated via the solution of the Nonlinear Schrödinger' ...
                ' Equation (NLSE) using the Split-Step Fourier Method (SSFM). Through ' ...
                'this simulation, the evolution of the pulse over distance is visualized, ' ...
                'offering insights into its temporal and spectral properties. By comparing ' ...
                'the parameters of the initial and final pulses, a deeper understanding of ' ...
                'how key factors such as fibre length, dispersion, and nonlinearity impact' ...
                ' the characteristics of the pulse as it traverses the fibre medium can ' ...
                'be attained.'],
                '',
                '',
                '\\begin{thebibliography}{9}',
                '\\bibitem{non_linear_optics}',
                ['``Nonlinear Fiber Optics'', Govind P.~Agrawal, 6$^{\\text{th}}$ edition,' ...
                ' Academic Press, New York (2019).'],
                '\\end{thebibliography}',
                '',
                '',
                '\\label{LastPage}',
                '\\end{document}',
            }];

            % compute number of lines in report cell array
            [num_rows,~] = size(report);
            %get set the file name 
            filename = 'EE4546_Assessment_1.tex'; 
            %save the file in a .tex file 
            file = fopen(filename, "wt");
            %write each line in the file
            for line_num = 1:num_rows
                fprintf(file, sprintf('%s\n',report{line_num}));
            end 
            fclose(file); 

            %display message to let the user know that a report has been
            %generated
            msgbox('EE4546_Assessment_1.tex report generated successfully.', 'Success', 'modal');

            %% Enable all buttons
            app.SimulateButton.Enable = "on";
            app.ReportButton.Enable = "on";
        end

        % Button pushed function: ExitButton
        function ExitButtonPushed(app, event)
            %close the app on exit button click
            delete(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create OpticalFibreSimulationUIFigure and hide until all components are created
            app.OpticalFibreSimulationUIFigure = uifigure('Visible', 'off');
            app.OpticalFibreSimulationUIFigure.Position = [100 100 786 549];
            app.OpticalFibreSimulationUIFigure.Name = 'OpticalFibreSimulation';

            % Create UIAxes
            app.UIAxes = uiaxes(app.OpticalFibreSimulationUIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Visible = 'off';
            app.UIAxes.Position = [350 10 418 282];

            % Create InputsLabel
            app.InputsLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.InputsLabel.FontSize = 16;
            app.InputsLabel.Interpreter = 'latex';
            app.InputsLabel.Position = [180 517 58 22];
            app.InputsLabel.Text = 'Inputs ';

            % Create OutputsLabel
            app.OutputsLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.OutputsLabel.FontSize = 16;
            app.OutputsLabel.Interpreter = 'latex';
            app.OutputsLabel.Position = [547 517 66 22];
            app.OutputsLabel.Text = 'Outputs';

            % Create DutycycleEditFieldLabel
            app.DutycycleEditFieldLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.DutycycleEditFieldLabel.HorizontalAlignment = 'right';
            app.DutycycleEditFieldLabel.Interpreter = 'latex';
            app.DutycycleEditFieldLabel.Position = [91 430 73 22];
            app.DutycycleEditFieldLabel.Text = 'Duty cycle:';

            % Create DutycycleEditField
            app.DutycycleEditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.DutycycleEditField.Limits = [1e-16 100];
            app.DutycycleEditField.AllowEmpty = 'on';
            app.DutycycleEditField.Position = [175 431 68 22];
            app.DutycycleEditField.Value = [];

            % Create Wavelengthlambda_0EditFieldLabel
            app.Wavelengthlambda_0EditFieldLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.Wavelengthlambda_0EditFieldLabel.HorizontalAlignment = 'right';
            app.Wavelengthlambda_0EditFieldLabel.Interpreter = 'latex';
            app.Wavelengthlambda_0EditFieldLabel.Position = [62 387 102 22];
            app.Wavelengthlambda_0EditFieldLabel.Text = 'Wavelength (\( \lambda_0 \)):';

            % Create Wavelengthlambda_0EditField
            app.Wavelengthlambda_0EditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.Wavelengthlambda_0EditField.Limits = [1e-16 1e+16];
            app.Wavelengthlambda_0EditField.AllowEmpty = 'on';
            app.Wavelengthlambda_0EditField.Position = [175 387 68 22];
            app.Wavelengthlambda_0EditField.Value = [];

            % Create DropDown_WavelengthUnit
            app.DropDown_WavelengthUnit = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_WavelengthUnit.Items = {'m', 'mm', 'μm', 'nm'};
            app.DropDown_WavelengthUnit.Position = [245 387 85 22];
            app.DropDown_WavelengthUnit.Value = 'μm';

            % Create FibrelengthLEditFieldLabel
            app.FibrelengthLEditFieldLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FibrelengthLEditFieldLabel.HorizontalAlignment = 'right';
            app.FibrelengthLEditFieldLabel.Interpreter = 'latex';
            app.FibrelengthLEditFieldLabel.Position = [63 339 101 22];
            app.FibrelengthLEditFieldLabel.Text = 'Fibre length (L):';

            % Create FibrelengthLEditField
            app.FibrelengthLEditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.FibrelengthLEditField.Limits = [1e-16 1e+16];
            app.FibrelengthLEditField.ValueDisplayFormat = '%11.5g';
            app.FibrelengthLEditField.AllowEmpty = 'on';
            app.FibrelengthLEditField.Position = [175 339 68 22];
            app.FibrelengthLEditField.Value = [];

            % Create DropDown_FibreLengthUnit
            app.DropDown_FibreLengthUnit = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_FibreLengthUnit.Items = {'m', 'km'};
            app.DropDown_FibreLengthUnit.Position = [245 339 85 22];
            app.DropDown_FibreLengthUnit.Value = 'km';

            % Create Lossalpha_textdBLabel
            app.Lossalpha_textdBLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.Lossalpha_textdBLabel.HorizontalAlignment = 'right';
            app.Lossalpha_textdBLabel.Interpreter = 'latex';
            app.Lossalpha_textdBLabel.Position = [94 299 70 22];
            app.Lossalpha_textdBLabel.Text = 'Loss (\( \alpha_{\text{dB}} \)):';

            % Create Lossalpha_textdBEditField
            app.Lossalpha_textdBEditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.Lossalpha_textdBEditField.Limits = [1e-16 1e+16];
            app.Lossalpha_textdBEditField.AllowEmpty = 'on';
            app.Lossalpha_textdBEditField.Position = [175 299 68 22];
            app.Lossalpha_textdBEditField.Value = [];

            % Create DropDown_LossUnit
            app.DropDown_LossUnit = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_LossUnit.Items = {'dB/m', 'dB/km'};
            app.DropDown_LossUnit.Position = [245 299 85 22];
            app.DropDown_LossUnit.Value = 'dB/km';

            % Create DispersioncoefficientDEditFieldLabel
            app.DispersioncoefficientDEditFieldLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.DispersioncoefficientDEditFieldLabel.HorizontalAlignment = 'right';
            app.DispersioncoefficientDEditFieldLabel.Interpreter = 'latex';
            app.DispersioncoefficientDEditFieldLabel.Position = [7 256 157 22];
            app.DispersioncoefficientDEditFieldLabel.Text = 'Dispersion coefficient (D):';

            % Create DispersioncoefficientDEditField
            app.DispersioncoefficientDEditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.DispersioncoefficientDEditField.Limits = [1e-16 1e+16];
            app.DispersioncoefficientDEditField.AllowEmpty = 'on';
            app.DispersioncoefficientDEditField.Position = [175 256 68 22];
            app.DispersioncoefficientDEditField.Value = [];

            % Create DropDown_DUnit
            app.DropDown_DUnit = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_DUnit.Items = {'s/m²', 'ps/nm/km'};
            app.DropDown_DUnit.Position = [245 256 85 22];
            app.DropDown_DUnit.Value = 'ps/nm/km';

            % Create Nonlinearcoefficientn_2Label
            app.Nonlinearcoefficientn_2Label = uilabel(app.OpticalFibreSimulationUIFigure);
            app.Nonlinearcoefficientn_2Label.HorizontalAlignment = 'right';
            app.Nonlinearcoefficientn_2Label.Interpreter = 'latex';
            app.Nonlinearcoefficientn_2Label.Position = [10 215 154 22];
            app.Nonlinearcoefficientn_2Label.Text = 'Nonlinear coefficient (\( n_2 \)):';

            % Create Nonlinearcoefficientn_2EditField
            app.Nonlinearcoefficientn_2EditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.Nonlinearcoefficientn_2EditField.Limits = [1e-16 1e+16];
            app.Nonlinearcoefficientn_2EditField.AllowEmpty = 'on';
            app.Nonlinearcoefficientn_2EditField.Position = [175 215 68 22];
            app.Nonlinearcoefficientn_2EditField.Value = [];

            % Create DropDown_NUnit
            app.DropDown_NUnit = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_NUnit.Items = {'m²/W', 'km²/W'};
            app.DropDown_NUnit.Position = [245 215 85 22];
            app.DropDown_NUnit.Value = 'km²/W';

            % Create EffectivecoreareaA_texteffLabel
            app.EffectivecoreareaA_texteffLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.EffectivecoreareaA_texteffLabel.HorizontalAlignment = 'right';
            app.EffectivecoreareaA_texteffLabel.Interpreter = 'latex';
            app.EffectivecoreareaA_texteffLabel.Position = [16 169 148 22];
            app.EffectivecoreareaA_texteffLabel.Text = 'Effective core area (\( A_{\text{eff}} \)):';

            % Create EffectivecoreareaA_texteffEditField
            app.EffectivecoreareaA_texteffEditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.EffectivecoreareaA_texteffEditField.Limits = [1e-16 1e+16];
            app.EffectivecoreareaA_texteffEditField.AllowEmpty = 'on';
            app.EffectivecoreareaA_texteffEditField.Position = [175 169 68 22];
            app.EffectivecoreareaA_texteffEditField.Value = [];

            % Create DropDown_AUnit
            app.DropDown_AUnit = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_AUnit.Items = {'m²', 'mm²', 'μm²'};
            app.DropDown_AUnit.Position = [245 169 85 22];
            app.DropDown_AUnit.Value = 'μm²';

            % Create SimulateButton
            app.SimulateButton = uibutton(app.OpticalFibreSimulationUIFigure, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);
            app.SimulateButton.Position = [109 89 100 23];
            app.SimulateButton.Text = 'Simulate';

            % Create ReportButton
            app.ReportButton = uibutton(app.OpticalFibreSimulationUIFigure, 'push');
            app.ReportButton.ButtonPushedFcn = createCallbackFcn(app, @ReportButtonPushed, true);
            app.ReportButton.Enable = 'off';
            app.ReportButton.Position = [238 89 100 23];
            app.ReportButton.Text = 'Report';

            % Create ExitButton
            app.ExitButton = uibutton(app.OpticalFibreSimulationUIFigure, 'push');
            app.ExitButton.ButtonPushedFcn = createCallbackFcn(app, @ExitButtonPushed, true);
            app.ExitButton.Position = [168 51 100 23];
            app.ExitButton.Text = 'Exit';

            % Create Label
            app.Label = uilabel(app.OpticalFibreSimulationUIFigure);
            app.Label.Position = [245 431 25 22];
            app.Label.Text = '%';

            % Create leftsqrtWrightLabel
            app.leftsqrtWrightLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.leftsqrtWrightLabel.Interpreter = 'latex';
            app.leftsqrtWrightLabel.Position = [674 476 45 22];
            app.leftsqrtWrightLabel.Text = '$\left[\sqrt{W}\right]$';

            % Create PeakamplitudeLabel
            app.PeakamplitudeLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PeakamplitudeLabel.Interpreter = 'latex';
            app.PeakamplitudeLabel.Position = [403 475 98 22];
            app.PeakamplitudeLabel.Text = 'Peak amplitude:';

            % Create PeakamplitudeLabel_val_initial
            app.PeakamplitudeLabel_val_initial = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PeakamplitudeLabel_val_initial.BackgroundColor = [1 1 1];
            app.PeakamplitudeLabel_val_initial.Position = [501 475 72 22];
            app.PeakamplitudeLabel_val_initial.Text = '';

            % Create psLabel_3
            app.psLabel_3 = uilabel(app.OpticalFibreSimulationUIFigure);
            app.psLabel_3.Interpreter = 'latex';
            app.psLabel_3.Position = [674 441 25 22];
            app.psLabel_3.Text = 'ps';

            % Create PulsepositionLabel
            app.PulsepositionLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PulsepositionLabel.Interpreter = 'latex';
            app.PulsepositionLabel.Position = [403 440 91 22];
            app.PulsepositionLabel.Text = 'Pulse position:';

            % Create PulsepositionLabel_val_initial
            app.PulsepositionLabel_val_initial = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PulsepositionLabel_val_initial.BackgroundColor = [1 1 1];
            app.PulsepositionLabel_val_initial.Position = [501 440 72 22];
            app.PulsepositionLabel_val_initial.Text = '';

            % Create psLabel_4
            app.psLabel_4 = uilabel(app.OpticalFibreSimulationUIFigure);
            app.psLabel_4.Interpreter = 'latex';
            app.psLabel_4.Position = [674 409 25 22];
            app.psLabel_4.Text = 'ps';

            % Create FWHMLabel
            app.FWHMLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FWHMLabel.Interpreter = 'latex';
            app.FWHMLabel.Position = [403 408 55 22];
            app.FWHMLabel.Text = 'FWHM:';

            % Create FWHMLabel_val_initial
            app.FWHMLabel_val_initial = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FWHMLabel_val_initial.BackgroundColor = [1 1 1];
            app.FWHMLabel_val_initial.Position = [501 408 72 22];
            app.FWHMLabel_val_initial.Text = '';

            % Create radps2Label_2
            app.radps2Label_2 = uilabel(app.OpticalFibreSimulationUIFigure);
            app.radps2Label_2.Interpreter = 'latex';
            app.radps2Label_2.Position = [674 374 47 22];
            app.radps2Label_2.Text = 'rad/ps$^2$';

            % Create ChirpLabel
            app.ChirpLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.ChirpLabel.Interpreter = 'latex';
            app.ChirpLabel.Position = [403 373 44 22];
            app.ChirpLabel.Text = 'Chirp:';

            % Create ChirpLabel_val_initial
            app.ChirpLabel_val_initial = uilabel(app.OpticalFibreSimulationUIFigure);
            app.ChirpLabel_val_initial.BackgroundColor = [1 1 1];
            app.ChirpLabel_val_initial.Position = [501 373 72 22];
            app.ChirpLabel_val_initial.Text = '';

            % Create radpsLabel_2
            app.radpsLabel_2 = uilabel(app.OpticalFibreSimulationUIFigure);
            app.radpsLabel_2.Interpreter = 'latex';
            app.radpsLabel_2.Position = [674 340 42 22];
            app.radpsLabel_2.Text = 'rad/ps';

            % Create FrequencyshiftLabel
            app.FrequencyshiftLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FrequencyshiftLabel.Interpreter = 'latex';
            app.FrequencyshiftLabel.Position = [403 339 99 22];
            app.FrequencyshiftLabel.Text = 'Frequency shift:';

            % Create FrequencyshiftLabel_val_initial
            app.FrequencyshiftLabel_val_initial = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FrequencyshiftLabel_val_initial.BackgroundColor = [1 1 1];
            app.FrequencyshiftLabel_val_initial.Position = [501 339 72 22];
            app.FrequencyshiftLabel_val_initial.Text = '';

            % Create radLabel_2
            app.radLabel_2 = uilabel(app.OpticalFibreSimulationUIFigure);
            app.radLabel_2.Interpreter = 'latex';
            app.radLabel_2.Position = [674 307 26 22];
            app.radLabel_2.Text = 'rad';

            % Create PhaseangleLabel
            app.PhaseangleLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PhaseangleLabel.Interpreter = 'latex';
            app.PhaseangleLabel.Position = [403 306 78 22];
            app.PhaseangleLabel.Text = 'Phase angle:';

            % Create PhaseangleLabel_val_initial
            app.PhaseangleLabel_val_initial = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PhaseangleLabel_val_initial.BackgroundColor = [1 1 1];
            app.PhaseangleLabel_val_initial.Position = [501 306 72 22];
            app.PhaseangleLabel_val_initial.Text = '';

            % Create PulsetypeDropDownLabel
            app.PulsetypeDropDownLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PulsetypeDropDownLabel.HorizontalAlignment = 'right';
            app.PulsetypeDropDownLabel.Interpreter = 'latex';
            app.PulsetypeDropDownLabel.Position = [123 129 70 22];
            app.PulsetypeDropDownLabel.Text = 'Pulse type:';

            % Create PulsetypeDropDown
            app.PulsetypeDropDown = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.PulsetypeDropDown.Items = {'Soliton', 'Gaussian'};
            app.PulsetypeDropDown.Position = [208 129 100 22];
            app.PulsetypeDropDown.Value = 'Soliton';

            % Create InitialLabel
            app.InitialLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.InitialLabel.Interpreter = 'latex';
            app.InitialLabel.Position = [520 496 41 22];
            app.InitialLabel.Text = 'Initial';

            % Create PeakamplitudeLabel_val_final
            app.PeakamplitudeLabel_val_final = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PeakamplitudeLabel_val_final.BackgroundColor = [1 1 1];
            app.PeakamplitudeLabel_val_final.Position = [590 476 72 22];
            app.PeakamplitudeLabel_val_final.Text = '';

            % Create PulsepositionLabel_val_final
            app.PulsepositionLabel_val_final = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PulsepositionLabel_val_final.BackgroundColor = [1 1 1];
            app.PulsepositionLabel_val_final.Position = [590 441 72 22];
            app.PulsepositionLabel_val_final.Text = '';

            % Create FWHMLabel_val_final
            app.FWHMLabel_val_final = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FWHMLabel_val_final.BackgroundColor = [1 1 1];
            app.FWHMLabel_val_final.Position = [590 409 72 22];
            app.FWHMLabel_val_final.Text = '';

            % Create ChirpLabel_val_final
            app.ChirpLabel_val_final = uilabel(app.OpticalFibreSimulationUIFigure);
            app.ChirpLabel_val_final.BackgroundColor = [1 1 1];
            app.ChirpLabel_val_final.Position = [590 374 72 22];
            app.ChirpLabel_val_final.Text = '';

            % Create FrequencyshiftLabel_val_final
            app.FrequencyshiftLabel_val_final = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FrequencyshiftLabel_val_final.BackgroundColor = [1 1 1];
            app.FrequencyshiftLabel_val_final.Position = [590 340 72 22];
            app.FrequencyshiftLabel_val_final.Text = '';

            % Create PhaseangleLabel_val_final
            app.PhaseangleLabel_val_final = uilabel(app.OpticalFibreSimulationUIFigure);
            app.PhaseangleLabel_val_final.BackgroundColor = [1 1 1];
            app.PhaseangleLabel_val_final.Position = [590 307 72 22];
            app.PhaseangleLabel_val_final.Text = '';

            % Create FinalLabel
            app.FinalLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.FinalLabel.Interpreter = 'latex';
            app.FinalLabel.Position = [609 497 37 22];
            app.FinalLabel.Text = 'Final';

            % Create TransmissionspeedEditFieldLabel
            app.TransmissionspeedEditFieldLabel = uilabel(app.OpticalFibreSimulationUIFigure);
            app.TransmissionspeedEditFieldLabel.HorizontalAlignment = 'right';
            app.TransmissionspeedEditFieldLabel.Interpreter = 'latex';
            app.TransmissionspeedEditFieldLabel.Position = [44 476 122 22];
            app.TransmissionspeedEditFieldLabel.Text = 'Transmission speed:';

            % Create TransmissionspeedEditField
            app.TransmissionspeedEditField = uieditfield(app.OpticalFibreSimulationUIFigure, 'numeric');
            app.TransmissionspeedEditField.Limits = [1e-16 1e+16];
            app.TransmissionspeedEditField.AllowEmpty = 'on';
            app.TransmissionspeedEditField.Position = [175 476 68 22];
            app.TransmissionspeedEditField.Value = [];

            % Create DropDown_TransmissionSpeed
            app.DropDown_TransmissionSpeed = uidropdown(app.OpticalFibreSimulationUIFigure);
            app.DropDown_TransmissionSpeed.Items = {'b/s', 'Kb/s', 'Mb/s', 'Gb/s'};
            app.DropDown_TransmissionSpeed.Position = [245 475 85 23];
            app.DropDown_TransmissionSpeed.Value = 'Gb/s';

            % Show the figure after all components are created
            app.OpticalFibreSimulationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FibreSimulationGUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.OpticalFibreSimulationUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.OpticalFibreSimulationUIFigure)
        end
    end
end