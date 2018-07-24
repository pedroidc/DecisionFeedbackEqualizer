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
M = 4;          % Ordem da modulação QPSK
NTbits = 1e8;   % Número total de bits a serem transmitidos
Nbits = 1e4;    % Número de bits a serem transmitidos por realização
Ntr = 1000;     % Comprimento da sequência de treinamento em símbolos
Ntrials = NTbits/Nbits; % Número total de realizações
Ns = Nbits/log2(M);
u = 0.01;

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
    eqm = zeros(Nbits, 1);
    Ph = 0;
    
    for iReal = 1:Ntrials
        
        %% Transmitter:
        xa = randsrc(Nbits, 1, [0 1]);          % Information bit stream
        xamod = step(hMod, xa);                 % Modulate
        
        %%  Channel
        hab = hba;                                % Reciprocal channel
        ybch = filter(hab, 1, xamod);             % Apply fading
        ybch = awgn(ybch, SNR(iSNR), 'measured'); % AWGN
        
        %% Receiver:
        % DFE:
        % Initializations:
        wf = zeros(Neq, 1); % Forward filter weights
        wb = zeros(Neq, 1); % Feedback filter weights
        yd = zeros(Ns, 1);  % DFE Output
        e = zeros(Ns, 1);   % DFE error
        yb = 0;             % Feedback filter output for first iteration
        for n = 1:Ns
            % Forward filter:
            uf = [ybch(n) uf(1:end-1)];    % Input vector
            yf =  uf*wf;                   % Output sample
            % Decision:
            yi = yf - yb;                           % Decision input sample
            yd(n) = step(hDemod, step(hMod, yi));   % Demodulated output
            % Feedback filter:
            ub = [yd(n) ub(1:end-1)];   % Input vector
            yb = ub*wb;                 % Output sample
            % Adaptation:
            e(n) = yi - ye;
            wf = wf(:, n) + mu*uf'*e(n);
            wb = wb(:, n) + mu*ub'*e(n);
        end
        
        % Demodulate:
        xbdemod = step(hDemod, yd);
        
        %% Performance Evaluation
        % BER in B:
        [~, bertemp] = biterr(xbdemod, xa);
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