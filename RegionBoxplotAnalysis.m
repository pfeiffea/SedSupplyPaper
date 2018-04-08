%% %%%%%%%%%%%
%
% The code runs the basic channel geometry and grain size analysis 
% described in Pfeiffer et al. (2017). 
%
% Files needed to run this code: 
% 'DataCompilation.txt' - compilation of channel geometry and grain data
% and
% 'Be10Data.mat' - compilation of 10Be erosion rates
%
% Allison Pfeiffer
%
%%%%%%%%%%%%%%

clear

%% Read in channel geometry and grain size data
id=fopen('DataCompilation.txt.');
M=textscan(id,'%n %s %*s %*s %*s %s %s %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',2,'Delimiter','\t'); 
fclose(id);

    S = M{5};
    D50 = M{6};
    h = M{7};
    w = M{8};
    DistPlate = M{12};

%% Read in 10Be data
id=fopen('10BeData.txt');
Be=textscan(id,'%*n %f %*s %f %*s %f','headerlines',1,'Delimiter','\t'); 
fclose(id);
 
Epub = Be{1}; % Original published erosion rate (mm/kyr)
E = Be{2}; % Willenbring's re-calculated erosion rate (mm/kyr)
EDistPlate = Be{3}; % Distance from the Pacific plate boundary or the west coast, whichever is further west... 
 
% % Willenbring didn't re-analyze Balco data (they were originally published using the preferred processing method)
% E(547:565) = Be10Data(547:565,2); %Ecalc = Epub


%% Determine region ('category') based on distance from plate boundary (calculated using Euclidian Distance tool in ArcGIS)

% Channel data
Category=ones(length(S),1); % Category = 1 
Category(DistPlate<=250000)=0; % Category = 0 if the sites are within 250km of the plate boundary

% Erosion rate data
ECategory=ones(length(E),1).*1;
ECategory(EDistPlate<=250000)=0;

%% Calculate Shields, etc
rho= 1000; % density of water
rhos=2650; % density of sediment
Rb= (rhos-rho)/rho; % bouyant density
g = 9.8; % gravitational acceleration
Rh = (w.*h)./(w+(2*h)); % estimated hydraulic radius

tauc = 0.15.*(S.^0.25); % Using Lamb et al. (2008) equation for estimating slope-dependent critical stress

% Calculate critical stress 
taubf = zeros(length(S),1);
for i=1:length(S)
    if isnan(w(i))== 1
        taubf(i)= (rho*g*h(i).*S(i))./((rhos-rho).*g.*D50(i));
    else
        taubf(i)= (rho*g*Rh(i).*S(i))./((rhos-rho).*g.*D50(i));
    end
end

%% Optional: testing alternative approaches to estimating tau*bf/tau*c

    % tauc = 2.18*S +0.021; % Mueller et al. (2005)
    % tauc = 0.39.*(S.^0.44)% Prancevic and Lamb (2015)
    % tauc = 0.36.*(S.^ 0.46); % Pitlick et al. (2008)
    % tauc = 0.035 % constant, not slope-dependent. 

    % Does hydraulic radius matter? Use depth instead  
    % taubf= (rho*g*h.*S)./((rhos-rho).*g.*D50); 

%% Ratio of bankfull Shields stress to critical Shields stress
Ratio = taubf./tauc;

%% Calculate Qt via Recking 2013 bedload transport model

% 0
D84 = 2.1.*D50;

% Eqn 12
tauSm = ((5.*S) + 0.06).*(D84./D50).^((4.4.*sqrt(S)) - 1.5);

% Eqn 13
for i=1:length(S)
    if isnan(w(i))== 1
        tauS84(i) = S(i).*h(i)./(1.65*D84(i)); 
    else
        tauS84(i) = S(i).*Rh(i)./(1.65*D84(i)); 
    end
end

% Eqn 11
Phi = 14.*(tauS84'.^2.5)./(1+(tauSm./tauS84').^4);

% Eqn ? 
qsv = sqrt(g.*1.65.*(D84.^3)).*Phi;

% 
Qt = 2650.*qsv;

 %% Boxplot Ratio -- Figure 1b in paper
figure(1) 

boxplot(Ratio,Category,'labels',{'West Coast','Other'},'colors','k','symbol','k.','outliersize',6)
% Test for statistically significant difference between region populations.
% H = 1--> we reject the null hypotheses that populations are statistically
% identical. p = p-value. 

[H,p] = ttest2(Ratio(Category==1),Ratio(Category==0),0.05,'both','unequal') % This is Welch's t-test. 

    % IMPROVEMENT not in 2017 manuscript: Welch's t-test should be used on data
    % with a normal distribution. Log normalize, then find p value. (Note: this
    % change results in a more significant difference between the populations).

    %[H,p] = ttest2(log10(Ratio(Category==1)),log10(Ratio(Category==0)),0.05,'both','unequal') % This is Welch's t-test. 

text(1.25,5,['p = ',num2str(roundn(p,-17))]);
text(.9,9,['n = ',num2str(length(Ratio(Category==0)))]);
text(.9,8,['median = ',num2str(roundn(median(Ratio(Category==0)),-3))]);
text(1.9,9,['n = ',num2str(length(Ratio(Category==1)))]);
text(1.9,8,['median = ',num2str(roundn(median(Ratio(Category==1)),-3))]);

ylabel('\tau*_{bf}/\tau*_{c}')
xlabel('Region')

%     
%% Erosion Boxplot -- Figure 1c in paper
figure(2)

boxplot(E,ECategory,'labels',{'West Coast','Other'},'colors','k','symbol','k.','outliersize',6)
% Test for statistically significant difference between region populations
[H,p] = ttest2(E(ECategory==0),E(ECategory==1),0.05,'both','unequal') % This is Welch's t-test. 

    % IMPROVEMENT not in 2017 manuscript: Welch's t-test should be used on data
    % with a normal distribution. Log normalize, then find p value. (Note: this
    % change results in a more significant difference between the populations).

    %[H,p] = ttest2(log10(E(ECategory==0)),log10(E(ECategory==1)),0.05,'both','unequal') % This is Welch's t-test. 

text(1.25,30,['p = ',num2str(roundn(p,-17))]);
text(.9,9,['n = ',num2str(length(E(ECategory==0)))]);
text(.9,6,['median = ',num2str(roundn(median(E(ECategory==0)),-3))]);
text(1.9,9,['n = ',num2str(length(E(ECategory==1)))]);
text(1.9,6,['median = ',num2str(roundn(median(E(ECategory==1)),-3))]);

    ylabel('Erosion rate (mm/kyr)')
    xlabel('Region')
    axis([.5 2.5 1 10000])
    set(gca,'yscale','log');
