%% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code analyzes the paired erosion and channel geometry sites in 
% Pfeiffer et al. (2017) to create Figure 2, therein. 
%
% Files needed to run this code:
% PairedSites.txt

% %%%%%%%%%%%%%%%%%%%%%%%%%

clear
clf

%% Load the data

id=fopen('PairedSites.txt');
M=textscan(id,'%*n %*s %*s %*s %*s %*s %f %f %f %*s %*s %*s %f','headerlines',2,'Delimiter','\t'); 
fclose(id);

Eshort = M{1};
E = M{2};
Ratio = M{3};
Region = M{4};

%% 10Be erosion rates
fig1 = figure(1);

x = log10(E);
y = log10(Ratio);
z = Region; 

% We only have short term erosion rates for some sites, ignore those. 
y(isnan(x)==1)=[];
z(isnan(x)==1)=[];
x(isnan(x)==1)=[];

plot(x(z==1),y(z==1),'k.', 'Color',[0 .4 .1],'MarkerSize',12); hold on
plot(x(z==0),y(z==0),'.', 'Color',[0.6 0.1 .2], 'MarkerSize',12);
ylabel('log(\tau*_{bf}/\tau*_c)')
xlabel('log(Erosion) log(mm/kyr)')

    % STATS
    coeffs=polyfit(x,y,1);
    % % Get fitted values
    fittedX = linspace(min(x), max(x), 10);
    fittedY = polyval(coeffs, fittedX);

    % % Plot the fitted line
    hold on;
    plot(fittedX, fittedY, 'k-', 'LineWidth', 1);

    [R,pval]=corrcoef(x,y);
    r2=R.^2;
    text(mean(x)*.5,max(y)*.9,['p= ',num2str(pval(2,1))]);
     text(mean(x)*.5,max(y)*.8,['r^2 ',num2str(r2(2,1))]);

%% Short term erosion rates

x = log10(Eshort);
y = log10(Ratio);
z = Region; 

% We only have 10Be erosion rates for some sites, ignore those. 
y(isnan(x)==1)=[];
z(isnan(x)==1)=[];
x(isnan(x)==1)=[];

plot(x(z==1),y(z==1),'ko','Color',[0 0.45 .1], 'MarkerSize',5); 
plot(x(z==0),y(z==0),'o','Color',[0.55 0.1 .2],'MarkerSize',5); hold on

    coeffs=polyfit(x,y,1);
    % % Get fitted values
    fittedX = linspace(min(x), max(x), 10);
    fittedY = polyval(coeffs, fittedX);

    % % Plot the fitted line
    hold on;
    plot(fittedX, fittedY, 'k--', 'LineWidth', 1);

    [R,pval]=corrcoef(x,y);
    r2=R.^2;
    text(mean(x)*.5,max(y),['pshort= ',num2str(pval(2,1))]);
     text(mean(x)*.5,max(y)*.95,['r^2 short ',num2str(r2(2,1))]);

%% Reference Lines showing median tau*bf/tau*c values, and median 10Be erosion rates for each region. 
% Note: I have hard-coded in the published medians. 

other = log10(1.03); hold on;
west = log10(2.35);
eother = log10(25);
ewest = log10(176);

oline = refline(0,other);
set(oline,'Color',[.7 .7 .7],'LineStyle','-','LineWidth', .75)
wline = refline(0,west);
set(wline,'Color',[.7 .7 .7],'LineStyle',':','LineWidth',.75 )

y1 = [-.4 3]; % arbitrary drawing limits for vertical reference lines. 
plot([eother eother], y1,'-','Color', [.7 .7 .7],'LineWidth',.75)
plot([ewest ewest], y1,':','Color', [.7 .7 .7],'LineWidth', .75)

ax = gca;
set(ax,'XTick',[.5 1 1.5 2 2.5 3])
set(ax,'YTick',[-.25 0 .25 .5 .75 1 ])
    
axis([.5 3 -.4 1])
axis square


hold off
