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
path7=[ '../TaskResrvoir/'];
addpath(genpath(path7));

load empirical_REST_dbs80.mat;
load hcp1003_REST1_LR_dbs80.mat;

N=62;
indexN=[1:31 50:80];
NSUB=1003;
Tmax=1200;

Isubdiag = find(tril(ones(N),-1));

TR=0.72;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

load SC_dbs80HARDIFULL.mat;
C = SC_dbs80HARDI;
C=C(indexN,indexN);
C=C-diag(diag(C));
C = C/max(max(C));

%%% FC empirical

for sub=1:NSUB
    sub
    clear signal_filt Phases;
    ts=subject{sub}.dbs80ts;
    ts=ts(indexN,:);
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
    end
    ts=signal_filt(:,50:end-50);
    FCe(sub,:,:)=corrcoef(ts');
end
FCemp=squeeze(mean(FCe));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear corrFC;
sigma=0.01;

p_inter=0.1;
p_interregion=0.1;
NQL=20;
kQL=10;
Nall=N*2*NQL;

Grange=0.01;%0.00005:0.00005:0.001;
Gloc=0.001;
% Gloc_range=0.02:0.02:0.2;

for i=1:N*NQL*2
    for b1=1:N
        if i<=b1*2*NQL && i>(b1-1)*2*NQL
            block(i)=b1;
        end
    end
end

for i=1:N*2*NQL
    f_diffall(i)=f_diff(block(i));
end
omega = 2*pi*f_diffall';

%%% Quantum-like

p_edge=0;
nnq=1;
for G=Grange

    wC=zeros(N*NQL*2,N*NQL*2);
    for i=1:N*NQL*2
        for j=1:N*NQL*2
            n1=block(i);
            n2=block(j);
            % if rand<p_interregion*C(n1,n2)
            %     wC(i,j)=G;
            % end
            wC(i,j)=G*C(n1,n2);
        end
    end
    for n=1:2*NQL:N*NQL*2
        CQL1=k_regular_graph(NQL,kQL,1000);
        CQL2=k_regular_graph(NQL,kQL,1000);
        for i=1:NQL
            for j=1:NQL
                if CQL1(i,j)==1 && rand<p_edge
                    CQL1(i,j)=0;
                end
                if CQL2(i,j)==1 && rand<p_edge
                    CQL2(i,j)=0;
                end
            end
        end
        CQL=zeros(NQL*2,NQL*2);
        CQL(1:NQL,1:NQL)=CQL1;
        CQL(1+NQL:end,1+NQL:end)=CQL2;
        for i=1:NQL
            for j=1+NQL:NQL*2
                if rand<p_inter
                    CQL(i,j)=1;
                    CQL(j,i)=1;
                end
            end
        end
        wC(n:NQL*2+n-1,n:NQL*2+n-1)=Gloc*CQL;
    end
    [FCsim2,COVsim,COVsimtotal,A]=hopf_int(wC,f_diffall,sigma);
    nn1=1;
    for n1=1:2*NQL:N*NQL*2
        nn2=1;
        for n2=1:2*NQL:N*NQL*2
            suma=[];
            for i=n1:n1+2*NQL-1
                for j=n2:n2+2*NQL-1
                    suma=[suma FCsim2(i,j)];
                end
            end
            FCsim(nn1,nn2)=mean(suma);
            nn2=nn2+1;
        end
        nn1=nn1+1;
    end

    %% fitting
    corrFC(nnq)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
    errFC(nnq)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)
    nnq=nnq+1;
end

%%% NON Quantum-like

p_edge=0.9;
nnq=1;
for G=Grange

    wC=zeros(N*NQL*2,N*NQL*2);
    for i=1:N*NQL*2
        for j=1:N*NQL*2
            n1=block(i);
            n2=block(j);
            % if rand<p_interregion*C(n1,n2)
            %     wC(i,j)=G;
            % end
            wC(i,j)=G*C(n1,n2);
        end
    end
    for n=1:2*NQL:N*NQL*2
        CQL1=k_regular_graph(NQL,kQL,1000);
        CQL2=k_regular_graph(NQL,kQL,1000);
        for i=1:NQL
            for j=1:NQL
                if CQL1(i,j)==1 && rand<p_edge
                    CQL1(i,j)=0;
                end
                if CQL2(i,j)==1 && rand<p_edge
                    CQL2(i,j)=0;
                end
            end
        end
        CQL=zeros(NQL*2,NQL*2);
        CQL(1:NQL,1:NQL)=CQL1;
        CQL(1+NQL:end,1+NQL:end)=CQL2;
        for i=1:NQL
            for j=1+NQL:NQL*2
                if rand<p_inter
                    CQL(i,j)=1;
                    CQL(j,i)=1;
                end
            end
        end
        wC(n:NQL*2+n-1,n:NQL*2+n-1)=Gloc*CQL;
    end
    [FCsim2,COVsim,COVsimtotal,A]=hopf_int(wC,f_diffall,sigma);
    nn1=1;
    for n1=1:2*NQL:N*NQL*2
        nn2=1;
        for n2=1:2*NQL:N*NQL*2
            suma=[];
            for i=n1:n1+2*NQL-1
                for j=n2:n2+2*NQL-1
                    suma=[suma FCsim2(i,j)];
                end
            end
            FCsim(nn1,nn2)=mean(suma);
            nn2=nn2+1;
        end
        nn1=nn1+1;
    end

    %% fitting
    corrFC(nnq)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
    errFC(nnq)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)
    nnq=nnq+1;
end

%%%%  Boxplots

%% QL

p_edge=0;

G=0.001;
for trials=1:100
    wC=zeros(N*NQL*2,N*NQL*2);
    for i=1:N*NQL*2
        for j=1:N*NQL*2
            n1=block(i);
            n2=block(j);
            if rand<p_interregion*C(n1,n2)
                wC(i,j)=G;
            end
            % wC(i,j)=G*C(n1,n2);
        end
    end
    for n=1:2*NQL:N*NQL*2
        CQL1=k_regular_graph(NQL,kQL,1000);
        CQL2=k_regular_graph(NQL,kQL,1000);
        for i=1:NQL
            for j=1:NQL
                if CQL1(i,j)==1 && rand<p_edge
                    CQL1(i,j)=0;
                end
                if CQL2(i,j)==1 && rand<p_edge
                    CQL2(i,j)=0;
                end
            end
        end
        CQL=zeros(NQL*2,NQL*2);
        CQL(1:NQL,1:NQL)=CQL1;
        CQL(1+NQL:end,1+NQL:end)=CQL2;
        for i=1:NQL
            for j=1+NQL:NQL*2
                if rand<p_inter
                    CQL(i,j)=1;
                    CQL(j,i)=1;
                end
            end
        end
        wC(n:NQL*2+n-1,n:NQL*2+n-1)=G*CQL;
    end
    [FCsim2,COVsim,COVsimtotal,A]=hopf_int(wC,f_diffall,sigma);
    nn1=1;
    for n1=1:2*NQL:N*NQL*2
        nn2=1;
        for n2=1:2*NQL:N*NQL*2
            suma=[];
            for i=n1:n1+2*NQL-1
                for j=n2:n2+2*NQL-1
                    suma=[suma FCsim2(i,j)];
                end
            end
            FCsim(nn1,nn2)=mean(suma);
            nn2=nn2+1;
        end
        nn1=nn1+1;
    end

    %% fitting
    corrFC_opt(trials)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
    errFC_opt(trials)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)
end

%% NON QL

p_edge=0.85;

G=0.0005;
for trials=1:100
    wC=zeros(N*NQL*2,N*NQL*2);
    for i=1:N*NQL*2
        for j=1:N*NQL*2
            n1=block(i);
            n2=block(j);
            if rand<p_interregion*C(n1,n2)
                wC(i,j)=G;
            end
            % wC(i,j)=G*C(n1,n2);
        end
    end
    for n=1:2*NQL:N*NQL*2
        CQL1=k_regular_graph(NQL,kQL,1000);
        CQL2=k_regular_graph(NQL,kQL,1000);
        for i=1:NQL
            for j=1:NQL
                if CQL1(i,j)==1 && rand<p_edge
                    CQL1(i,j)=0;
                end
                if CQL2(i,j)==1 && rand<p_edge
                    CQL2(i,j)=0;
                end
            end
        end
        CQL=zeros(NQL*2,NQL*2);
        CQL(1:NQL,1:NQL)=CQL1;
        CQL(1+NQL:end,1+NQL:end)=CQL2;
        for i=1:NQL
            for j=1+NQL:NQL*2
                if rand<p_inter
                    CQL(i,j)=1;
                    CQL(j,i)=1;
                end
            end
        end
        wC(n:NQL*2+n-1,n:NQL*2+n-1)=G*CQL;
    end
    [FCsim2,COVsim,COVsimtotal,A]=hopf_int(wC,f_diffall,sigma);
    nn1=1;
    for n1=1:2*NQL:N*NQL*2
        nn2=1;
        for n2=1:2*NQL:N*NQL*2
            suma=[];
            for i=n1:n1+2*NQL-1
                for j=n2:n2+2*NQL-1
                    suma=[suma FCsim2(i,j)];
                end
            end
            FCsim(nn1,nn2)=mean(suma);
            nn2=nn2+1;
        end
        nn1=nn1+1;
    end

    %% fitting
    corrFC_opt_non(trials)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
    errFC_opt_non(trials)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)
end
save results_HopfQLbit.mat corrFC corrFC_non corrFC_opt corrFC_opt_non;