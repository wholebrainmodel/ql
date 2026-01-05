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
dtt=0.1;

clear corrFC errFC;
sigma=0.01;
sig=sigma;
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

p_interregion=0.1;
NQL=40;
kQL=20;
Nall=N*NQL;

D = 0.5*(sigma^2)*eye(2*Nall);
iD = inv(D);

Grange=0.01:0.01:0.2;
Gloc_range=1;%0.0005:0.0005:0.01;

for i=1:N*NQL
    for b1=1:N
        if i<=b1*NQL && i>(b1-1)*NQL
            block(i)=b1;
        end
    end
end

for i=1:N*NQL
    f_diffall(i)=f_diff(block(i));
end
omega = 2*pi*f_diffall';

p_edge=0.2; %% good -> 0, bad-> 0.95

nnq=1;
for Gloc=Gloc_range
    mmq=1;
    for G=Grange
        wC=zeros(N*NQL,N*NQL);
        for i=1:N*NQL
            for j=1:N*NQL
                n1=block(i);
                n2=block(j);
                if rand<p_interregion*C(n1,n2)
                    wC(i,j)=G;
                    wC(j,i)=G;
                end
            end
        end
        for n=1:NQL:N*NQL
            CQL=k_regular_graph(NQL,kQL,1000);
            for i=1:NQL
                for j=1:NQL
                    if CQL(i,j)==1 && rand<p_edge
                        CQL(i,j)=0;
                    end
                end
            end
            wC(n:NQL+n-1,n:NQL+n-1)=Gloc*CQL;
        end
        [FCsim2,COVsim,COVsimtotal,A]=hopf_int(wC,f_diffall,sigma);
        nn1=1;
        for n1=1:NQL:N*NQL
            nn2=1;
            for n2=1:NQL:N*NQL
                suma=[];
                for i=n1:n1+NQL-1
                    for j=n2:n2+NQL-1
                        suma=[suma FCsim2(i,j)];
                    end
                end
                FCsim(nn1,nn2)=mean(suma);
                nn2=nn2+1;
            end
            nn1=nn1+1;
        end

        %% fitting
        corrFC(nnq,mmq)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
        errFC(nnq,mmq)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)
        mmq=mmq+1;
    end
    nnq=nnq+1;
end

%%%%  Boxplots

%% QL

p_edge=0.2;

[aux ind]=min(squeeze(errFC(1,:)));
G=Grange(ind);

for trials=1:100
    wC=zeros(N*NQL,N*NQL);
    for i=1:N*NQL
        for j=1:N*NQL
            n1=block(i);
            n2=block(j);
            if rand<p_interregion*C(n1,n2)
                wC(i,j)=G;
                wC(j,i)=G;
            end
        end
    end
    for n=1:NQL:N*NQL
        CQL=k_regular_graph(NQL,kQL,1000);
        for i=1:NQL
            for j=1:NQL
                if CQL(i,j)==1 && rand<p_edge
                    CQL(i,j)=0;
                end
            end
        end
        wC(n:NQL+n-1,n:NQL+n-1)=Gloc*CQL;
    end

    %%
        [V DD]=eig(wC);
        DD=diag(real(DD));
        [dmax]=sort(DD,'descend');
        numql=length(find(dmax>3*std(DD)));
        spectralgap_ql(trials)=sum(abs(diff(dmax(1:numql+1))));
    %%

    [FCsim2,COVsim,Theta,J]=hopf_int(wC,f_diffall,sigma);
    nn1=1;
    for n1=1:NQL:N*NQL
        nn2=1;
        for n2=1:NQL:N*NQL
            suma=[];
            for i=n1:n1+NQL-1
                for j=n2:n2+NQL-1
                    suma=[suma FCsim2(i,j)];
                end
            end
            FCsim(nn1,nn2)=mean(suma);
            nn2=nn2+1;
        end
        nn1=nn1+1;
    end
    %% fitting
    corrFC_ql(trials)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
    errFC_ql(trials)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)

    %% Energy
    % A=-J;
    % EntroFlow_ql(trials)=abs(trace(A'*iD*A*Theta-A))
    %%
    A=-J;
    Theta=(Theta+Theta')/2;
    Theta=Theta+dt*(-(A*Theta+Theta*A')+2*D);
    DTA=D*inv(Theta)-A;
    [P L]=eig(Theta);
    AP=A*P;
    TT=P'*Theta*P;
    XTX=AP*TT*AP';
    for node=1:2*Nall
        EntroFlow_ql_nodes_tr(trials,node)=iD(node,node)*XTX(node,node)-A(node,node);
    end
    EntroFlow_ql(trials)=mean(abs(squeeze(EntroFlow_ql_nodes_tr(trials,:))));
    %%
end

EntroFlow_ql_nodes2=squeeze(mean(EntroFlow_ql_nodes_tr));
nn1=1;
for n1=1:NQL:N*NQL
    suma=[];
    for i=n1:n1+NQL-1
        suma=[suma EntroFlow_ql_nodes2(i)];
    end
    EntroFlow_ql_nodes(nn1)=mean(suma);
    nn1=nn1+1;
end

%% NON QL
p_edge=0.8; %% good -> 0, bad-> 0.95

for trials=1:100
    wC=zeros(N*NQL,N*NQL);
    for i=1:N*NQL
        for j=1:N*NQL
            n1=block(i);
            n2=block(j);
            if rand<p_interregion*C(n1,n2)
                wC(i,j)=G;
                wC(j,i)=G;
            end
        end
    end
    for n=1:NQL:N*NQL
        CQL=k_regular_graph(NQL,kQL,1000);
        for i=1:NQL
            for j=1:NQL
                if CQL(i,j)==1 && rand<p_edge
                    CQL(i,j)=0;
                end
            end
        end
        wC(n:NQL+n-1,n:NQL+n-1)=Gloc*CQL;
    end

    %%
    [V DD]=eig(wC);
    DD=diag(real(DD));
    [dmax]=sort(DD,'descend');
    numql=length(find(dmax>3*std(DD)));
    spectralgap_non(trials)=sum(abs(diff(dmax(1:numql+1))));
    
    %%

    [FCsim2,COVsim,Theta,J]=hopf_int(wC,f_diffall,sigma);
    nn1=1;
    for n1=1:NQL:N*NQL
        nn2=1;
        for n2=1:NQL:N*NQL
            suma=[];
            for i=n1:n1+NQL-1
                for j=n2:n2+NQL-1
                    suma=[suma FCsim2(i,j)];
                end
            end
            FCsim(nn1,nn2)=mean(suma);
            nn2=nn2+1;
        end
        nn1=nn1+1;
    end
    %% fitting
    corrFC_non(trials)=corr2(FCemp(Isubdiag),FCsim(Isubdiag))
    errFC_non(trials)=mean((FCemp(Isubdiag)-FCsim(Isubdiag)).^2)

    %% Energy
    % A=-J;
    % EntroFlow_non(trials)=abs(trace(A'*iD*A*Theta-A))
        %%
    A=-J;
    Theta=(Theta+Theta')/2;
    Theta=Theta+dt*(-(A*Theta+Theta*A')+2*D);
    DTA=D*inv(Theta)-A;
    [P L]=eig(Theta);
    AP=A*P;
    TT=P'*Theta*P;
    XTX=AP*TT*AP';
    for node=1:2*Nall
        EntroFlow_non_nodes_tr(trials,node)=iD(node,node)*XTX(node,node)-A(node,node);
    end
    EntroFlow_non(trials)=mean(abs(squeeze(EntroFlow_non_nodes_tr(trials,:))));
    %%
end
EntroFlow_non_nodes2=squeeze(mean(EntroFlow_non_nodes_tr));
nn1=1;
for n1=1:NQL:N*NQL
    suma=[];
    for i=n1:n1+NQL-1
        suma=[suma EntroFlow_non_nodes2(i)];
    end
    EntroFlow_non_nodes(nn1)=mean(suma);
    nn1=nn1+1;
end

ranksum(corrFC_non,corrFC_ql)
ranksum(errFC_non,errFC_ql)
ranksum(spectralgap_non,spectralgap_ql)
ranksum(EntroFlow_non,EntroFlow_ql)

figure(1)
violinplot([corrFC_non' corrFC_ql'])
figure(2)
violinplot([errFC_non' errFC_ql'])
figure(3)
violinplot([spectralgap_non' spectralgap_ql'])
figure(4)
violinplot([EntroFlow_non' EntroFlow_ql'])

save results_HopfQLinterregion4020_pedge0208.mat corrFC_ql corrFC_non errFC_ql errFC_non spectralgap_non spectralgap_ql ...
    EntroFlow_non EntroFlow_ql EntroFlow_non_nodes EntroFlow_ql_nodes G;