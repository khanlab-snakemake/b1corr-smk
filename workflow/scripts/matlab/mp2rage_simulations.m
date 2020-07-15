clc;
clear all;

%% MP2RAGE parameters for simulations
MP2RAGE.London.B0=7;           % in Tesla
MP2RAGE.London.TR=6;           % MP2RAGE TR in seconds
MP2RAGE.London.TRFLASH=6.7e-3; % TR of the GRE readout
MP2RAGE.London.TIs=[800e-3 2700e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.London.NZslices=[56 112];% Slices Per Slab * [PartialFourierInSlice-0.5  0.5]
MP2RAGE.London.FlipDegrees=[4 5];% Flip angle of the two readouts in degrees
MP2RAGE.London.PartialFourierInSlice=0.75

MP2RAGE.Maastricht.B0=7;
MP2RAGE.Maastricht.TR=5;
MP2RAGE.Maastricht.TRFLASH=6.9e-3;
MP2RAGE.Maastricht.TIs=[900e-3 2750e-3];
MP2RAGE.Maastricht.NZslices=[120 120];
MP2RAGE.Maastricht.FlipDegrees=[5 3];
MP2RAGE.Maastricht.PartialFourierInSlice=1

MP2RAGE.ProtocolA.B0=7;
MP2RAGE.ProtocolA.TR=6;
MP2RAGE.ProtocolA.TRFLASH=6.9e-3;
MP2RAGE.ProtocolA.TIs=[800e-3 2700e-3];
MP2RAGE.ProtocolA.NZslices=[56 112];
MP2RAGE.ProtocolA.FlipDegrees=[4 5];
MP2RAGE.ProtocolA.PartialFourierInSlice=0.75

%% Define signal and noise functions as in the original MP2RAGE paper
% (https://doi.org/10.1016/j.neuroimage.2009.10.002)
Signalres = @(x1,x2) x1.*x2./((x2.^2+x1.^2));
noiseres = @(x1,x2) ((x2.^2-x1.^2).^2 ./(x2.^2+x1.^2).^3 ).^(0.5);

%% Typical T1 values at 7T
T1WM=1.1;
T1GM=1.85;
T1CSF=3.9;

%% T1 and B1 ranges for simulations
T1range=0:.05:T1CSF;
B1range=0:.1:2;
FlipAngles=combvec(1:1:15,1:1:15);

%% Simulated data
MP2RAGEamps=zeros(size(B1range,2),100,numel(fieldnames(MP2RAGE)),size(FlipAngles,2));
T1vectors=zeros(size(B1range,2),100,numel(fieldnames(MP2RAGE)),size(FlipAngles,2));
CNRs=zeros(size(B1range,2),size(T1range,2),numel(fieldnames(MP2RAGE)),size(FlipAngles,2));

k=0;
for site={'London','Maastricht','Maastricht-2'}
    k=k+1;
    PARAMS=MP2RAGE.(site{1})
    
    l=0;
    for B1=B1range
        l=l+1;
        
        n=0;
        for FlipAngle=FlipAngles
            n=n+1;
            
            [MP2RAGEamp T1vector IntensityBeforeComb]=MP2RAGE_lookuptable(2,PARAMS.TR,PARAMS.TIs,B1*FlipAngle',PARAMS.NZslices,PARAMS.TRFLASH,'normal',[],1);
            MP2RAGEamps(l,:,k,n)=MP2RAGEamp;
            T1vectors(l,:,k,n)=T1vector;
            
            m=0;
            for T1=T1range
                m=m+1;
                
                [temp, preT1 ]=min(abs(T1 - T1vector));
                [temp, postT1 ]=min(abs((T1+0.05) - T1vector));
                
                Signal=Signalres(IntensityBeforeComb([preT1,postT1],1),IntensityBeforeComb([preT1,postT1],2));
                Noise=noiseres(IntensityBeforeComb([preT1,postT1],1),IntensityBeforeComb([preT1,postT1],2));
                CNRs_PFcorr(l,m,k,n) = 1000 * sum(((Signal(2:end)-Signal(1:(end-1)))./sqrt(Noise(2:end).^2+Noise(1:(end-1)).^2))./sqrt(PARAMS.TR)).*sqrt(PARAMS.PartialFourierInSlice);
            end
        end
    end
end

save('/home/ROBARTS/rhaast/graham/scratch/B1corr/output/simulations/data_PFcorrected.mat','MP2RAGEamps','T1vectors','CNRs_PFcorr','T1range','B1range','FlipAngles');

%% Plot B1+ sensitivity
k=0;
for B1=B1range
    k=k+1;
    [MP2RAGEamp T1vector IntensityBeforeComb]=MP2RAGE_lookuptable(2,MP2RAGE.ProtocolA.TR,MP2RAGE.ProtocolA.TIs,B1*MP2RAGE.ProtocolA.FlipDegrees,MP2RAGE.ProtocolA.NZslices,MP2RAGE.ProtocolA.TRFLASH,'normal',[],1);
    plot(MP2RAGEamp,T1vector,'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
end