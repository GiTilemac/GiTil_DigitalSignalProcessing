%% Frequency Estimation using Eigenanalysis of the Autocorrelation Matrix
%
% Author: Tilemachos S. Doganis
%
%% Description
% This script simulates N realizations of a stochastic process of the form
% S_omega(n) = A*e^(j*omega*n + phi), n = 0,1,...,M.
%
% omega: The process is of rank P, meaning it contains P separate
% frequencies 'omega'
% A: Signal amplitudes, A_p = 1/(2^p), p = 0,1,...,P-1
% phi: Phase offset, uniformly distributed in [-pi,pi], mutually
%      independent between realizations. This parameter represents the
%      random start time of each signal in the process.
%
%   The purpose of the script is to estimate signal and noise parameters,
% mainly the signal frequencies, using the Pisarenko Harmonic Decomposition,
% MUltiple SIgnal Classification (MUSIC) and EigenVector (EV) algorithms.
%   A precondition for the above methods is that the rank of the stochastic
% process, i.e. the number of the desired frequencies, is known a priori.
%   The estimation process is repeated for different numbers of realizations,
% as well as different sized submatrices of the stochastic process
% autocorrelation matrix 'Rxx'.
%
%% Process steps
% 1. Produce N realizations of the stochastic process, with duration M,
% resulting in a MxN matrix 'S'. Add AWGN 'W' with power sigma_w = 0.5,
% so that X = S + W.
%
% 2. Calculate the autocorrelation matrix 'R_XX' of 'X', and its eigenvalues
% 'lambda' and eigenvectors 'U' (MxM).
%
% Since AWGN is stationary in the wide sense, and independent from the
% signal, its autocorrelation matrix 'R_WW' will be diagonal and additive
% to the autocorrelation matrix of the signal 'R_SS'. Given that the noise
% mean is zero, 'R_WW' will contain the signal power in the diagonal.
% Therefore, it follows that: R_XX = |A|^2*e_M*e_M^h + sigma_W^2,
% where e_M(w) = [ 1 e^(-jw) e^(-2jw) ... e^-((M-1)jw) ]
% By multiplying the autocorrelation matrix with vector e_M, it follows:
% R_XX*e_M = (M*|A|^2 + sigma_w^2)*e_M, meaning that e_M is an eigenvector
% of R_XX, corresponding to the eigenvalue M*|A|^2 + sigma_w^2.
% Considering the noise power is low enough, the first 'P' eigenvectors
% span the signal space 'U_S', while the other 'M-P' the noise space
% 'U_N'. The same observation applies to the corresponding eigenvalues.
% Since the basis vectors for these spaces are orthogonal to each other,
% <e_M,n_m> = 0, m = P+1,...,M. Thus, R_XX*n_m = sigma_W^2*n_m =
% lambda_m*n_m.
%
% 3. Estimate the noise variance using the arithmetic mean of the last 'M-P'
% eigenvalues, which correspond to the noise space.
%
% 4. Define the trigonometric polynomials described in the Pisarenko
% Harmonic Decomposition method: P_(M,m)(e^(jw)) = e_M(w)^h*n_m,
% m = P+1, P+2, ..., M, where n_m are the noise eigenvectors. Then define
% the rational functions: Q_(M,m)(e^jw) = 1/|P_(M,m)(e^jw)|^2
% The Q function is infinite at the frequencies 'w' where
% <e_M(w),n_m> = 0, meaning that these signals do not belong the noise
% space, and consequently are a part of the signal space.
% 
% 5. Find the most common roots of the polynomials P_(M,m), using MATLAB's
% roots() function and k-means clustering. The biggest P clusters contain
% the root estimations as centers.
%
% 6. By multiplying the autocorrelation matrix left and right with the
% the signal space vectors and solving the ensuing equation for 'A', the
% signal amplitudes can be estimated.
%
% 7. Define the polynomials according to the Pisarenko Harmonic Decomposition,
% MUSIC and EV algorithms in order to plot them. The estimated frequencies are
% detected using the top peaks.
%
% The above process is repeated multiple times for different parameters
% 'M' and 'N', in order to study the effectiveness of the three algorithms.
% According to the results of this script, MUSIC seems to be the most
% effective one in this case.
%
%% Program code
clear
clc
close all

% Define stochastic process & other parameters
N_reals = [ 100 1000 10000 50000];   % # of realizations
M = [10 20 30 40 50];                % Max. Rank of Autocorrelation
P = 5;                               % Stochastic Process rank
sigma_w = 0.5;                       % Noise power
A = 1./(2.^(0:4)');                  % Signal Amplitudes
omega = [0.2 0.4 0.5 0.75 0.88]'*pi; % Signal frequencies (angular)

%% Initialization
N = N_reals(end);                    % Include all realizations at begin.
S = zeros(M(end),N);                 % Initialize Signal vectors
W = sqrt(sigma_w).*randn(M(end),N);  % Noise vectors
n = repmat((0:M(end)-1)',1,N);       % Discrete time vectors
phi = -pi + 2*pi*rand(P,N);          % Phase offsets

% Produce and sum up P stochastic signals S_i = A_i*e^(j*w_i*n + phi_i)
% (M time units x N realizations)
for j = 1:P
    S = S + A(j).*( exp(1i*( omega(j).*n + repmat(phi(j,:),[M(end) 1]) )));
end

% Add AWGN
X_full = S + W;

%% Estimation
for Nind = 1:4
    
    % Keep N_i realizations
    N = N_reals(Nind); 
    X = X_full(:,1:N);
    fprintf(['Running Eigenanalysis with ',num2str(N),' realizations...\n']);
    
    % Parameter estimation matrix R_M (for P = 5 separate processes)
    R_M = cell(P,1);
    for k = 1:P
        R_M{k} = struct('e_M',[],'A',zeros(P,1),'sigma',[],'omega',zeros(P,1),'lambda',...
            [],'cten1',[],'cten2',[],'P_Mm',[],'Q_Mm',[],'Q_MUSIC',[],'kw',[]);
    end

    % Estimate stochastic mean
    mX = 1/N*sum(X,2);

    % Estimate autocorrelation matrix
    Rxx = X*X';
    
    % Repeat process for Rxx submatrices (reverse order to compute Rxx once)
    for k = 5:-1:1;
        % Eigenanalysis on submatrix of Rxx
        Rxx = Rxx(1:M(k),1:M(k));    % Slice submatrix of Rxx
        [~,lambda,U] = svd(Rxx);     % Apply SVD
        lambda = diag(lambda)/N;     % Normalize eigenvalues
        R_M{k}.lambda = lambda;
        
        % Normalize signal eigenvectors (P first eigenvectors), so that the
        % real element on the first place is the unit.
        U(:,1:P) = bsxfun( @rdivide,U(:,1:P),U(1,1:P)); 
        
        % Signal space is spanned by first P eigenvectors, noise space
        % is spanned by the M - P remaining ones
        U_S = U(:,1:P); 
        U_N = U(:,(P+1):end); 
        
        % Noise eigenvalues
        lambda_w = lambda((P+1):end);    
        
        %% 2.1) Estimate noise variance
        R_M{k}.sigma = mean(lambda_w);
        
        %% 2.3) 1st & 2nd order central tendencies
        R_M{k}.cten1 = 1/M(k)*sum(   lambda_w  -  mean(lambda_w) );     % Estimate 1st-order central tendency
        R_M{k}.cten2 = 1/M(k)*sum(  (lambda_w  -  mean(lambda_w) ).^2); % Estimate 2nd-order central tendency
        
        %% 4) Define trigonometric polynomials P_(M,m)(e^jw)
        R_M{k}.P_Mm = U_N; % Polynomial coefficients are noise eigenvectors
        
        %% 5) Estimate w_i's using common roots of all P_(M,m)
        rts = zeros(M(k)-1,M(k)-P);  % rts: roots of P_M,m(e^jw)

        % Estimate roots of P_(M,m) & solve for frequency 'w'
        for i = 1:(M(k)-P)
            rts(:,i) =  real(-1i*log(roots(R_M{k}.P_Mm(:,i))));    
        end
        rts = rts(:);                % Store roots in 1 vector
        rts = rts(rts>0 & rts<=pi);  % Keep only roots in [0,pi]
        n_roots = M(k)-P;            % #roots = #noise eigenvectors (M-P)

        % Apply k-means on all roots where #cluster = #roots
        [idx, kw, sumd] = kmeans(rts(:),n_roots);      
        
        %% Common Root Selection Algorithm
         
        % Find the clusters with the P largest populations
        cluster_pop = zeros(n_roots,1);
        for i = 1:n_roots
            cluster_pop(i) = sum(idx==i);
        end
        max_pop = sort(cluster_pop,'descend');
        max_pop = max_pop(1:P);
        
        % Keep all common roots of all polynomials, if none, keep roots
        % who solve the most
        criterium = min(min(max_pop),n_roots);  
        k_omegas = kw(cluster_pop >= criterium);
        
        % If more than P potential frequencies remain, re-cluster to P
        % clusters to unite neighboring centers
        if numel(k_omegas) > P 
            [~, k_omegas] = kmeans(k_omegas,P);   
        end
        
        % The cluster centers represent the frequency estimations
        R_M{k}.omega = sort(k_omegas,'ascend');
        
        %% 3) Estimate Signal Amplitudes
        % Bring the P frequencies w_P back into the form e^(j*w_P*n)
        e_P = exp(-1i.*R_M{k}.omega*(0:M(k)-1))'; 
        
        p = diag(U_S'*Rxx*U_S)/N;
        Mu = abs(U_S'*e_P).^2;
        R_M{k}.A = real(sqrt(Mu\(p-R_M{k}.sigma)));
        R_M{k}.A = sort(R_M{k}.A,'descend');
        
        %% 6) Define Rational Functions Q_(M,m)(e^jw)
        R_M{k}.kw = sort(kw);
        e_M = exp(-1i*R_M{k}.kw*(0:M(k)-1))';   % Evaluate e_M at cluster means, kw
        P_m = e_M'*R_M{k}.P_Mm;                 % Elements of n_m are coefficients of P_Mm
        R_M{k}.Q_Mm = 1./abs(P_m).^2;           % Define Rational Functions Q_(M,m)

        %% 7) Define Q_MUSIC
        R_M{k}.Q_MUSIC = 1./sum(abs(P_m).^2,2);

        %% 8) Define Q_EV 
        R_M{k}.Q_EV = 1./sum( bsxfun(@times, abs(P_m).^2,  (1./lambda_w)) ,2);
        
    end
    figure
    for k = 1:5
        subplot(1,5,k),histogram(R_M{k}.lambda(2:end),'BinMethod','sturges'),title(['\lambda_m(M=', num2str(M(k)), ')']);
    end
    figure,
    for k = 1:5
    %% 6) Plot Rational Functions Q_(M,m)(e^jw)
     subplot(3,5,k),hold on
    for i=1:(M(k)-P)
        plot(R_M{k}.kw,R_M{k}.Q_Mm(:,i)),xlim([0 3.2]), 
    end
    title(['Q_{(',num2str(M(k)),',m)}(e^{j\omega)}']),xlabel('\omega');
    hold off

    %% 7) Plot Q_MUSIC
    subplot(3,5,k+5),plot(R_M{k}.kw,R_M{k}.Q_MUSIC),xlim([0 3.2]),
    title(['Q_{',num2str(M(k)),'}^{MUSIC}(e^{j\omega)}']),xlabel('\omega');

    %% 8) Plot Q_EV
    subplot(3,5,k+10),plot(R_M{k}.kw,R_M{k}.Q_EV),xlim([0 3.2]),
    title(['Q_{',num2str(M(k)),'}^{EV}(e^{j\omega)}']),xlabel('\omega');
    end
    
    % Print results
    print_subm(R_M,A,sigma_w,omega);

end