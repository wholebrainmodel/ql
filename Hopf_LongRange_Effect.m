clear all;
path2=[ '../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../Tenet/TENET/'];
addpath(genpath(path3));
path4=[ '../TaskResrvoir/DataHCP100ordered'];
addpath(genpath(path4));
path5=[ '../Turbulence/Basics'];
addpath(genpath(path5));
path6=[ '../Turbulence/SC_longrange'];
addpath(genpath(path6));

load (['hcp971finalorder_REST1_LR_schaefer116.mat']);

N=100;
indexN=1:N;
sigma=0.01;
NSUB=971;
Isubdiag = find(tril(ones(N),-1));

%% FCemp
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

%% Meta empirical

for nsub=1:NSUB
    nsub
    ts=subject{nsub}.schaeferts;  % fMRI
    ts=ts(indexN,:);
    clear signal_filt Phases;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    Phases=Phases(:,50:end-50);
    Metaemp(nsub)=std(abs(sum(exp(1i*Phases),1))/N);
end

%% Spectral gaps.

load results_Ceff_Rest_sub116.mat;

Ceff=Ceff_rest_sub;
for nsub=1:NSUB
    [V D]=eig(squeeze(Ceff(nsub,:,:)));
    D=diag(real(D));
    [dmax]=sort(D,'descend');
    numql=length(find(dmax>3*std(D)));
    spectralgap(nsub)=sum(abs(diff(dmax(1:numql+1))));
end

[cc pp]=corr(Metaemp',spectralgap')
scatter(Metaemp',spectralgap','filled')

%%% get rid of Longrange
load COG_Schaefer100.mat;
dd=cog;

for i=1:N
    for j=1:N
        distancia(i,j)=sqrt(sum((dd(i,:)-dd(j,:)).^2));
        if distancia(i,j)>80
            dmask(i,j)=0;
        else
            dmask(i,j)=1;
        end
    end
end

for nsub=1:NSUB
    Ceffmod=squeeze(Ceff(nsub,:,:));
    Ceffmod=Ceffmod.*dmask;
    [V D]=eig(Ceffmod);
    D=diag(real(D));
    [dmax]=sort(D,'descend');
    numql=length(find(dmax>3*std(D)));
    spectralgapmod(nsub)=sum(abs(diff(dmax(1:numql+1))));
end

[cc pp]=corr(Metaemp',spectralgapmod')
scatter(Metaemp',spectralgapmod','filled')

ranksum(spectralgap,spectralgapmod)
boxplot([spectralgap' spectralgapmod']);

save results_spectralgap_meta_lonrange.mat spectralgap spectralgapmod Metaemp;

