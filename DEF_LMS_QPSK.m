%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of a communication system using decision-feedback equalizer
% (DFE) with the LMS algorithm in order to mitigate the frequency selective
% effects of the wireless channel.
%
% Author: Pedro Ivo da Cruz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; %close all;

%% Simulation Parameters:
M = 4;                  % QPSK modulation order
NTbits = 1e6;           % Total number bits to be transmitted
Nbits = 1e4;            % Number of bits for each trial
Ntr = 1000;             % Length of the training sequence
Ntrials = NTbits/Nbits; % Number of trials necessary
Ns = Nbits/log2(M);     % Number os symbols in each trial
mu = 0.01;              % LMS adaptation step
Neq = 13;               % DFE filters lenght

%% QPSK Modulation:
hMod = comm.RectangularQAMModulator( ...
    'ModulationOrder', M, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
hDemod = comm.RectangularQAMDemodulator( ...
    'ModulationOrder', M, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

%% Channels parameters:
SNR = 0:2:30;           % Signal-to-noise Ratio for simulation
L = 6;                  % Channel length
pp = exp(-(0:L-1)/2);   % Channel power profile with exponential decay with factor 2
pp = pp.'/sum(pp);      % Normalize to make the total power 1;

%% Initializations:
hError = comm.ErrorRate('ResetInputPort', true);
NSNR = length(SNR);
berB = zeros(NSNR, 1);
berE = zeros(NSNR, 1);
ecma = zeros(Ns, 1);

%% Simulation:
cont = 0;
barStr = 'SIM_QPSK_THP.m';
progbar = waitbar(0, 'Initializing waitbar...', ...
    'Name', barStr);

fprintf('Simulação iniciada\n')

for iSNR = 1:NSNR
    
    berBtrial = 0;
    berEtrial = 0;
    eqm = zeros(Ns, 1);
    Ph = 0;
    
    for iReal = 1:Ntrials
        
        %% Transmitter:
        x = randsrc(Nbits, 1, [0 1]); % Information bit stream
        xmod = step(hMod, x);         % Modulate
        
        %%  Channel
        h = sqrt(0.5*pp).*(randn(L, 1) + 1j*randn(L, 1)); % Fading Channel
        ych = filter(h, 1, xmod);                         % Apply fading
        ych = awgn(ych, SNR(iSNR), 'measured');           % AWGN
        
        %% Receiver:
        % DFE:
        % Initializations:
        wf = zeros(Neq, 1); % Forward filter weights
        wb = zeros(Neq, 1); % Feedback filter weights
        yd = zeros(Ns, 1);  % DFE Output
        e = zeros(Ns, 1);   % DFE error
        yb = 0;             % Feedback filter output for first iteration
        uf = zeros(1, Neq); % Forward filter input vector
        ub = zeros(1, Neq); % Feedback filter input vector
        for n = 1:Ns
            % Forward filter:
            uf = [ych(n) uf(1:end-1)]; % Input vector
            yf =  uf*wf;               % Output sample
            % Decision:
            yi = yf - yb;                           % Decision input sample
            yd(n) = step(hMod, step(hDemod, yi));   % Demodulated output
            % Feedback filter:
            ub = [yd(n) ub(1:end-1)]; % Input vector
            yb = ub*wb;               % Output sample
            % Adaptation:
            e(n) = yi - yb;
            wf = wf + mu*uf'*e(n);
            wb = wb + mu*ub'*e(n);
        end
        
        % Demodulate:
        hDemod.release();
        xbdemod = step(hDemod, yd);
        
        %% Performance Evaluation
        % BER in B:
        [~, bertemp] = biterr(xbdemod, x);
        berB(iSNR) = berB(iSNR) + bertemp;
        
        % DFE convergence:
        eqm = eqm + e.^2;
        
        % Simulation progress:
        cont = cont + 1;
        % Update bar:
        prog = cont/(Ntrials*NSNR);
        perc = 100*prog;
        waitbar(prog, progbar, sprintf('%.2f%% Concluído', perc));
        
    end
    
end

berB = berB/Ntrials;
eqm = eqm/Ntrials;

close(progbar)

EbNo = SNR - 10*log10(log2(M));
% BER teórica de um sistema QAM
% berTheory = berawgn(EbNo, 'qam', M);
berTheory = berfading(EbNo, 'qam', M, 1);
% berTheory = berawgn(EbNo, 'psk', M, 'nondiff');

% Plota os resultados
figure
semilogy(SNR, berB, 'bo--');
hold on;
semilogy(SNR, berE, 'ro--');
semilogy(SNR, berTheory, 'kx-');
legend('Alice', 'Eve', 'Theory');
xlabel('SNR (dB)');
ylabel('BER');
%axis([EbNo(1) EbNo(end) 0.0001 1]);
grid on;

figure
plot(10*log10(eqm));

% scatterplot(yach)

% fileStr = ['Resultados/THPrecoding_EveMODAdder_Nch' num2str(L)];

% Salva dados para posterior utilização:
% save(fileStr, ...
%     'berB', 'berE', 'ecma', 'SNR', 'EbNo', 'L');