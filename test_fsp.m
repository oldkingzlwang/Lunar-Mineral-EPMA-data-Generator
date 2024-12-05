function StrctFrm=test_fsp(data)
% Modified from Jesse B. Walters's program: MinPlot (Mineralogia, vol. 53, 2022, 51-66)

[~,n]=size(data);
for k = 1:n
    data(isnan(data(:, k)), k) = 0;
end
MC=zeros(size(data));
SiO2=60.083;TiO2=79.865;Al2O3=101.961;
Cr2O3=151.989;FeO=71.8442;MnO=70.937;
MgO=40.304;CaO=56.0774;Na2O=61.979;
K2O=94.195;P2O5=141.943;
W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O,P2O5];
for i=1:size(W,2)
    if i==3 || i==4 || i==9 || i==10 || i==11
        MC(:,i)=data(:,i)./W(:,i).*2;
    else
        MC(:,i)=data(:,i)./W(:,i);
    end
end
Opfu=8.0; %oxygens per formula unit

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Cr2O3
O2(:,5)=MC(:,5); %for FeO
O2(:,6)=MC(:,6); %for MnO
O2(:,7)=MC(:,7); %for MgO
O2(:,8)=MC(:,8); %for CaO
O2(:,9)=MC(:,9)./2; %for Na2O
O2(:,10)=MC(:,10)./2; %for K2O
O2(:,11)=MC(:,11).*(5/2); %for P2O5

O2total=sum(O2,2); %O2 totals
MCnormfact=Opfu./O2total; %normalization factor

%% Structural Formula
StrctFrm=MCnormfact.*MC; %creates a matrix of normalized cations
StrctFrm(:,12)=sum(StrctFrm,2); %calculations the total, which should be close to 5
