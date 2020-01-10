%% Script for CodeBook10 by Connor Haines 

    %%Summary 
    
    %This codebook uses the Black-Litterman model to create a buy-and-hold
    %portfolio that aims to maximize sharpe ratio. In order to form views
    %inputed into the Black-Litterman model, I first use Arma-Garch forecasting. 
    % Using these forecasts and additional discretionary views, I form  views that I
    %then input into the Black-Litterman model. I use the output from
    %Black-Litterman to maximize sharpe ratio. 
    
    %%Step-By-Step Overview 
    
    %This script uses imported price data from money.net to create a buy and
    %hold, fully invested, long-only investment portfolio made up of US
    %mid-cap stocks.
    
    %I used discretion to cull my investable universe to 20 possible
    %stocks. These are all stocks that feature prominintely in IJH, the
    %benchmark portfolio. 
    
    %As per my investor identity, my portfolio is optimized based on Sharpe
    %ratio.
    
    %My portfolio is formed at the beginning of the investmen horizon
    %(this semester) and held until the present. 
    
    %The Import Data section of this script imports
    %data from money.net.
    
    %The Portfolio Construction section uses the Black-Litterman model to 
    %construct my portfolio. Views are formed using a combination of
    %Arma-Garch forecasting and discretion.
    
    %The performance measurement section measures several performance
    %metrics. This section is adapted from code by Gonzalo Asis (Jan. 2017).
    
    %The performance report section constructs a ppt visually displaying my
    %portfolio performance. My results are compared to the benchmark
    %portfolio and a backwards-looking max sharpe ratio portfolio. 

%% Housekeeping
    clear all % clear all variables from memory
    close all % close all figures
    clc % clear command window

%% User Inputs

    %Data dates 
    startDate={'12-Jan-2017'}; %start of investment horizon 
    endDate={'19-Apr-2017'}; %end of investment horizon 
    startAnalysis={'12-Jan-2007'}; %start of backtesting period 
    endAnalysis={'11-Jan-2017'}; %end of backtesting period 
    
    %PowerPoint Slides 
    PresName = 'HainesConnor_PortfolioSummary10.pptx'; % name of ppt file
    reportNumber = 'Performance Report 10'; %number of codebook
    macroCommentary = {'Housing starts fell sharply this month (1.215m vs 1.262m consensus), as did manufacturing production (-0.4% M/M), leading to further doubt over the strength of the economy during the first quarter.', '10-yr treasury yields have fallen sharply over the last week to 2.18%.', 'Retail sales fell 0.2% in March, while the February number was revised four tenths lower.', 'CPI fell 0.3% in March, while core CPI fell 0.1%.'}; %macro commentary for this week 
    marketCommentary = {'UK Prime Minister May has announced she aims to hold the next UK  election on June 8th, well ahead of the next scheduled election in 2020. Analysts say she is hoping for a convincing win to strengthen the conservative’s position in parliament. The House of Commons will vote on the proposition today.', 'The Dow posted a sharpe 0.6% loss to 20,404 on Tuesday, but stocks were mixed with the Nasdaw up 0.2% to 5,863.', 'Oil fell by around 4% ($2 a barrel) on Wednesday, it’s most in 6-weeks.', 'Earnings season kicks off this week, with generally high expectations for a profitable quarter.'}; %market commentary for this week
    
 %% Import Data 
 %Creates two strucures. One called 'historical' containting historical
 %closing prices for all securities considered over backtesting period
 %and one called 'holdingperiod' with closing prices for all securites
 %considered over holding period. These could be pulled all at once, but it
 %is easier to use them later when they are stored in separate structures. 


 
 %Money.net Connection
 % Input username and password
    username = 'unc7@unc.edu';
    pwd = 'moneynet';
    
    c = moneynet(username,pwd); % create connection
    
    
 % Retrieve historical prices for securities
    symbols = {'RJF','SNPS','RMD','INGR', 'AMD', 'SIVB','ATO', 'PKG', 'ARE', 'DRE', 'Y', 'SBNY', 'DPZ', 'IT', 'STLD', 'RE', 'ANSS', 'CDNS', 'ALGN', 'UGI'};
    date = [datetime(startAnalysis) datetime(endAnalysis)];
    interval = '1D'; % time interval
    f = {'Close'}; % data fields we want
    % Need to loop through individual retrieval
    for j = 1:length(symbols) % loop through the length of symbols cells
        symbol = symbols{j}; % select symbol
        historical.(char(symbol)) = timeseries(c,symbol,date,interval,f); % retrieve data
    end
    
 % Retrieve holding period prices for  securities
    symbols = {'RJF','SNPS','RMD','INGR', 'AMD', 'SIVB','ATO', 'PKG', 'ARE', 'DRE', 'Y', 'SBNY', 'DPZ', 'IT', 'STLD', 'RE', 'ANSS', 'CDNS', 'ALGN', 'UGI'};
    date = [datetime(startDate) datetime(endDate)];
    interval = '1D'; % time interval
    f = {'Close'}; % data fields we want
    % Need to loop through individual retrieval
    for j = 1:length(symbols) % loop through the length of symbols cells
        symbol = symbols{j}; % select symbol
        holdingperiod.(char(symbol)) = timeseries(c,symbol,date,interval,f); % retrieve data
    end 
    
    % Risk-free rate  
        rf = 0.000164111; % daily risk-free from 1-month t-bill 2/8/17
        riskfreevariance=0; %sort of an assumption since t-bils aren't technically risk free
    
    % Retrieve historical prices for benchmark
    %this could be done with original datapull, but it's easier to navigate
    %when it's kept separate 
    date = [datetime(startAnalysis) datetime(endAnalysis)];
    interval = '1D'; % time interval
    f = {'Close'}; % data fields we want
    mrkthistorical = timeseries(c,'IJH',date,interval,f); % retrieve data
    
    %Retrieve holding period prices for benchmark
    date = [datetime(startDate) datetime(endDate)];
    mrktholding = timeseries(c,'IJH',date,interval,f); % retrieve data

        
    close(c); % close connection
    
 
%% Daily Returns for each Security and Benchmark   
    
    clearvars -except historical holdingperiod Reflection macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber rf riskfreevariance symbols

     % Compute Daily Returns for each Security and Benchmark
    
    for i=1:length(symbols)   
        symbol=symbols{i};
        Ret(:,i) = historical.(char(symbol)){2:end,2}./historical.(char(symbol)){1:end-1,2}-1; %compute daily returns
    end
    bmr(:,:)=mrkthistorical{2:end,2}./mrkthistorical{1:end-1,2}-1; %compute daily returns of benchmark

   
%% Forecast using ARMA GARCH to establish views for my portfolio
 
   %Size the data
    [T,N] = size(Ret); 

%Estimate during the estimation window
    %Set the model assuming an ARMA(1,1) GARCH(1,1) structure
        Mdl = arima('ARLags',1,'MALags',1,'Variance',garch(1,1)); 
    %Set the forecast horizon
        ForHor = 31; 
    %Loop through each asset
    forecastret =[];
    for i = 1:N
        clear ret 
    %Form the temporary returns
        ret = Ret(1:end,i); 
    %Estimate the model 
        Options=optimoptions(@fmincon,'ConstraintTolerance',1e-1,'Algorithm','sqp');
        EstMdl = estimate(Mdl,ret,'Options',Options);
    %Infer the variances
        [E0,V0,LogL] = infer(EstMdl,ret); 
    %Forecast
        [rfor,YMSE,V] = forecast(EstMdl,ForHor,'Y0',ret,'E0',E0,'V0',V0);
    %Append that forecast to the Return Series
        forecastret(:,i) = rfor;
       
    %Display where you are in the loop 
        clc
        td = ['Asset #',num2str(i)];
        disp(td); 
    end %End loop through assets
            
    
%% Set Views (Using Forecast) 

%This creates my pick matrix. I'm setting an absolute view for each asset
%based on its expected mean return from my forecasts, so I just need a
%diagonal matrix of ones. Note, you can also set relative views, but all my
%views are absolute. 
views=ones(N,1);
views=diag(views);

%This is where I set the expected returns for each view (the 'Q' matrix in
%the BL specifications). I get the expected returns from my forecast. 
viewexpectations=mean(forecastret)'; 

%Additional Discretionary Views 
viewexpectations(5,1)=viewexpectations(5,1).*2; %discretionary view that AMD will return above forecast
viewexpectations(18,1)=viewexpectations(18,1).*2; %discretionary view that CDNS will return above forecast
viewexpectations(19,1)=viewexpectations(19,1).*2; %discretionary view that ALGN will return above forecast


 %% Portfolio Construction using Black Litterman 
 
%BL settings

%risk aversion calibrated outside model.  Value can be changed for
%robustness
delta = 3; 

%Equilibrium weights in market portfolio.
weq = []; 
weq(1,1:N)=1/N;  
%The exact weights have changed slighty throughout the semester, but as each stock in my discretionary universe is held fairly prominently in IJH it works as an approximation to have each included equally as the equilibrium configuration. 
weq = weq+(1-sum(weq))/size(weq,1); %ensure weq sums to 1. 


sigma = cov(Ret);%prior covariance matrix
tau = 1/(T-N); %uncertainty in prior estimate of mean returns; common calibration. 
P = views; %"Pick" matrix
Q = viewexpectations;
Omega = diag(diag(P*tau*sigma*P')); %Variances of views; Common to scale with prior variances
[PosteriorReturn,PosteriorCov, PosteriorW, PosteriorPW] = hlblacklitterman(delta, weq, sigma, tau, P, Q, Omega); 

%Maximize SR and grab the weights
p = Portfolio('AssetMean',PosteriorReturn, 'AssetCovar',PosteriorCov); 
p = setDefaultConstraints(p);
swgt = estimateMaxSharpeRatio(p);
%In Sample risk and return estimates
[srsk,sret] = estimatePortMoments(p,swgt); 
    
clearvars -except benchmark SampleStats startingCash Reflection student_name weights Blotter bmr riskfreevariance rf swgt srsk sret hwgt hrsk hret cret crsk mret mrsk p pret prsk pwgt Ret historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols

    
%% Backwards Looking Max Sharpe Ratio Portfolio (for comparisson)
%Forms an sharpe-ratio optimized portfolio by backtesting the last ten
%years of historical prices. Utilizes the Financial Tools toolbox. 
    
% Compute Mean Return for each Security 
    AssetMean=mean(Ret)'; %average of all the daily returns for each security 
    bmmean=mean(bmr); %average daily return for benchmark 
    

% Compute Return Covariances 
    AssetCovar=cov(Ret);
    BmVar=var(bmr);
   
%Set variables for Portfolio Element     
mret = bmmean; %mean return of benchmark
mrsk = sqrt(BmVar); %std. dev. of bechmark returns
cret = rf; %risk free rate
crsk = riskfreevariance; %risk free=assume no variance 
    
% Create Portfolio Element     
    p = Portfolio('AssetList', symbols, 'RiskFreeRate', cret);
    p = setAssetMoments(p, AssetMean, AssetCovar);
    p = setDefaultConstraints(p); %requires fully-invested long-only portfolio 
    pwgt = estimateFrontier(p, 20); %construct efficient frontier for graph later
   [prsk, pret] = estimatePortMoments(p, pwgt); %construct efficient frontier for graph later

% Create Maximize the Sharpe Ratio Portfolio
% A maximum sharpe ratio portfolio seeks to maximize return per unit of
% risk

p = setInitPort(p, 0);

hwgt = estimateMaxSharpeRatio(p); %stores weight allocation for each security
[hrsk, hret] = estimatePortMoments(p, hwgt); %returns avg. daily risk and return of maximum sharpe ratio portfolio over previous two years 

clearvars -except benchmark SampleStats startingCash Reflection student_name weights swgt srsk sret hwgt hrsk hret cret crsk mret mrsk p Ret historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols

%% Performance Measurement 
%This section records a variety of performance measurement metrics for my
%portfolio. Holding period is 1/12/17 to
%4/19/17. 

%% Inputs
student_name = 'Connor Haines';
weights = num2cell(swgt(:,1:end)); %convert weights array to cell
benchmark='IJH';
startingCash = 1000000;

%% Performance of individual securities

% Create structure to hold sample statistics
if i==1
SampleStats = struct;
end

% Compute holding period returns for each security 
for i=1:length(symbols)
    symbol=symbols{i};
    SampleStats.ret_hold(i) = holdingperiod.(char(symbol)){end,2}./holdingperiod.(char(symbol)){1,2}-1;
end

% Compute Daily Returns for each Security 
    ret=[];
    for i=1:length(symbols)   
        symbol=symbols{i};
        ret(:,i) = holdingperiod.(char(symbol)){2:end,2}./holdingperiod.(char(symbol)){1:end-1,2}-1; %compute daily returns
    end
    SampleStats.ret = ret;



% Compute sample moments for the returns
SampleStats.median_ret = median(SampleStats.ret);
SampleStats.mean_ret = mean(SampleStats.ret);
voltemp=std(SampleStats.ret);
SampleStats.vol_ret = voltemp;
skewtemp=skewness(SampleStats.ret);
SampleStats.skew_ret = skewtemp;
kurtosistemp = kurtosis(SampleStats.ret);
SampleStats.kurt_ret = kurtosistemp;

% Compute cumulative returns over the period
cumsumtemp=cumsum(ret);
SampleStats.cumret=cumsumtemp;

clearvars -except benchmark SampleStats startingCash Reflection student_name weights Blotter swgt srsk sret hwgt hrsk hret cret crsk mret mrsk p pret prsk pwgt Ret historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols


%% Performance of Benchmark
% Create array of benchmark portfolio values
benchmarkshares=startingCash/mrktholding{1,2};
benchmarkport=benchmarkshares.*mrktholding{1:end,2};

% Compute returns of benchmark 
ret_mkt=mrktholding{2:end,2}./mrktholding{1:end-1,2}-1; %compute daily returns of benchmark
cumret_mkt=cumsum(ret_mkt); % cumulative returns

% Add to SampleStats structure
SampleStats.ret_mkt = ret_mkt;
SampleStats.cumret_mkt = cumret_mkt;

% Compute sample moments for benchmark returns
SampleStats.median_benchmarkport = median(ret_mkt);
SampleStats.mean_benchmarkport = mean(ret_mkt);
SampleStats.vol_benchmarkport = std(ret_mkt);
SampleStats.skew_benchmarkport = skewness(ret_mkt);
SampleStats.kurt_benchmarkport = kurtosis(ret_mkt);

% Sharpe ratio benchmark
SampleStats.benchmarksharpe = sharpe(ret_mkt,cret);


% Maximum drawdown benchamrk
SampleStats.benchmarkmax_draw = maxdrawdown(benchmarkport);

% Minimum and Maximum values of benchmark
SampleStats.benchmarkmax = max(benchmarkport);
SampleStats.benchmarkmin = min(benchmarkport);

clearvars -except ret_mkt benchmarkshares benchmarkport Reflection benchmark SampleStats startingCash student_name weights Blotter hwgt hret hrsk swgt srsk sret cret crsk mret mrsk p pret prsk pwgt Ret historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols

%% Portfolio Performance - loop through BL portfolio and backwards looking max sharpe portfolio 


%index of portfolios 
portindex={swgt; hwgt};

for i=1:length(portindex)

% Determine how many shares of each to buy
shares=[];
for j=1:length(symbols)
    symbol=symbols{j};
    shares(1,j) =portindex{i}(j,1).*startingCash./holdingperiod.(char(symbol)){1,2};
end

% Construct portfolio
for j=1:length(symbols)
    symbol=symbols{j};
    sec(:,j) = shares(j).*holdingperiod.(char(symbol)){:,2};
end

port(:,i)=sum(sec,2);

end


%performance metrics for portfolios
% Daily returns
SampleStats.ret_bl = price2ret(port(:,1));
SampleStats.ret_ms = price2ret(port(:,2));

% Cumulative returns
SampleStats.cumret_bl = cumsum(SampleStats.ret_bl);
SampleStats.cumret_ms = cumsum(SampleStats.ret_ms);

% Compute sample moments for the returns
SampleStats.median_bl = median(SampleStats.ret_bl);
SampleStats.mean_bl = mean(SampleStats.ret_bl);
SampleStats.vol_bl = std(SampleStats.ret_bl);
SampleStats.skew_bl = skewness(SampleStats.ret_bl);
SampleStats.kurt_bl = kurtosis(SampleStats.ret_bl);
SampleStats.median_ms = median(SampleStats.ret_ms);
SampleStats.mean_ms = mean(SampleStats.ret_ms);
SampleStats.vol_ms = std(SampleStats.ret_ms);
SampleStats.skew_ms = skewness(SampleStats.ret_ms);
SampleStats.kurt_ms = kurtosis(SampleStats.ret_ms);

% Volatility
SampleStats.volbl = std(port(:,1)); % total volatility
SampleStats.volbl5 = std(port(end-4:end,1)); % 5-day volatility
SampleStats.volms = std(port(:,2)); % total volatility
SampleStats.volms5 = std(port(end-4:end,2)); % 5-day volatility


% Sharpe ratio
SampleStats.sharpe_bl = sharpe(SampleStats.ret_bl, cret);
SampleStats.sharpe_ms = sharpe(SampleStats.ret_ms, cret); 

% Maximum drawdown
SampleStats.maxdraw_bl = maxdrawdown(port(:,1));
SampleStats.maxdraw_ms = maxdrawdown(port(:,2));

% Minimum and Maximum values
SampleStats.max_bl = max(port(:,1));
SampleStats.min_bl = min(port(:,1));
SampleStats.max_ms = max(port(:,2));
SampleStats.min_ms = min(port(:,2));

clearvars -except symbol ret_mkt ret_port benchmarkshares Reflection benchmarkport benchmark SampleStats startingCash student_name weights hwgt hret hrsk swgt srsk sret cret crsk mret mrsk p pret prsk pwgt Ret historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols


%% Additoinal Risk Metrics 

%Compute the Historical VaR for each portfolio
%Also create and save figure for ppt 
confidence_level = 0.95;
plot_flag = 1;
figure
subplot(3,1,1) 
SampleStats.histVaR_bl= -computeHistoricalVaR(SampleStats.ret_bl,confidence_level,plot_flag);  
subplot(3,1,2)
SampleStats.histVaR_benchmarkport= -computeHistoricalVaR(SampleStats.ret_mkt,confidence_level,plot_flag); 
subplot(3,1,3)
SampleStats.histVaR_ms= -computeHistoricalVaR(SampleStats.ret_ms,confidence_level,plot_flag); 
saveas(gcf,'Hist_var.png');

%Compute the Gaussian VaR,lpm for each 
RiskThreshold = 1-confidence_level;
PortValue = 1; 
SampleStats.gaussianVaR_bl = portvrisk(SampleStats.mean_bl,SampleStats.vol_bl,RiskThreshold,PortValue); 
SampleStats.lowerPartial_bl= lpm(SampleStats.ret_bl); 
SampleStats.gaussianVaR_ms = portvrisk(SampleStats.mean_ms,SampleStats.vol_ms,RiskThreshold,PortValue); 
SampleStats.lowerPartial_ms= lpm(SampleStats.ret_ms);
SampleStats.gaussianVaR_benchmarkport = portvrisk(SampleStats.mean_benchmarkport,SampleStats.vol_benchmarkport,RiskThreshold,PortValue); 
SampleStats.lowerPartial_benchmarkport= lpm(SampleStats.ret_mkt);


%% Dates
date = datenum(holdingperiod.(char(symbol)){1:end,1});

%% Create Figures
% Pie chart of asset allocation - black litterman
for i=1:length(symbols)
    labels(i) = strcat(symbols(i),' (',num2str(swgt(i)*100),'%)');
end

% Asset allocation pie chart - black litterman 
figure
pie((swgt),labels)
title('Asset Allocation - Black Litterman')
saveas(gcf,'Pie_alloc_bl.png');

% Histogram of returns - number of days return of portfolio falls into each
% 'bucket'. 
figure
histogram(SampleStats.ret_bl)
title('Returns Histogram')
xlabel('Daily Return in % buckets')
ylabel('Number of Days')
saveas(gcf,'Hist_port.png');

%compare to benchmark

% Daily return of portfolio and benchmark
%update: graph fixed so that returns start at zero
figure
date = [date(1,1)-1; date];
retporttemp=SampleStats.ret_bl;
retporttemp=[0; retporttemp];
plot(date(1:end-1,:),retporttemp) % shorten date because we lose 1 data point when computing yield
hold on
retmkttemp = SampleStats.ret_mkt;
retmkttemp = [0; retmkttemp];
plot(date(1:end-1,:),retmkttemp)
benchmark='IJH';
legend('Portfolio',char(benchmark))
title('Daily Yield')
ylabel('Return in %')
datetick('x',23,'keeplimits')
saveas(gcf,'Daily_ret.png');

% Cumulative return of portfolio and benchmark
figure
cumretporttemp = SampleStats.cumret_bl;
cumretporttemp = [0; cumretporttemp];
plot(date(1:end-1,:),cumretporttemp)
hold on
cumretmkttemp = SampleStats.cumret_mkt;
cumretmkttemp = [0; cumretmkttemp];
plot(date(1:end-1,:),cumretmkttemp)
legend('Portfolio',char(benchmark))
title('Cumulative Returns')
ylabel('Return in %')
datetick('x',23,'keeplimits') % place date in x-axis in date format
saveas(gcf,'Cum_ret.png');

%compare to backwards looking max sharpe portfolio

% Pie chart of asset allocation - backwards max sharpe
for i=1:length(symbols)
    labels2(i) = strcat(symbols(i),' (',num2str(hwgt(i)*100),'%)');
end

% Asset allocation pie chart - black litterman 
figure
subplot(2,1,1)
pie((swgt),labels)
title('Black Litterman')
subplot(2,1,2)
pie((hwgt),labels2)
title('Backwards Max Sharpe')
saveas(gcf,'Pie_alloc_ms.png');

% Daily return of portfolio and benchmark
%update: graph fixed so that returns start at zero
figure
retporttemp=SampleStats.ret_bl;
retporttemp=[0; retporttemp];
plot(date(1:end-1,:),retporttemp) % shorten date because we lose 1 data point when computing yield
hold on
retmstemp = SampleStats.ret_ms;
retmstemp = [0; retmstemp];
plot(date(1:end-1,:),retmstemp)
benchmark='Backwards Max Sharpe';
legend('Portfolio',char(benchmark))
title('Daily Yield')
ylabel('Return in %')
datetick('x',23,'keeplimits')
saveas(gcf,'Daily_ret_ms.png');

% Cumulative return of portfolio and benchmark
figure
cumretporttemp = SampleStats.cumret_bl;
cumretporttemp = [0; cumretporttemp];
plot(date(1:end-1,:),cumretporttemp)
hold on
cumretmstemp = SampleStats.cumret_ms;
cumretmstemp = [0; cumretmstemp];
plot(date(1:end-1,:),cumretmstemp)
legend('Black Litterman Portfolio',char(benchmark))
title('Cumulative Returns')
ylabel('Return in %')
datetick('x',23,'keeplimits') % place date in x-axis in date format
saveas(gcf,'Cum_ret_ms.png');
   %% Portfolio Construction   
    
    clearvars -except date benchmarkshares benchmarkport Reflection benchmark SampleStats startingCash student_name weights hwgt hret hrsk swgt srsk sret cret crsk mret mrsk p pret prsk pwgt Ret historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols


%% Create Tables
numbers = zeros(5,2);

%BL compard to Benchmark
% Fill out the table with the data from SampleStats
numbers(1,1)=SampleStats.cumret_bl(end,:);
numbers(1,2)=SampleStats.cumret_mkt(end,1);
numbers(2,1)=SampleStats.sharpe_bl;
numbers(2,2)=SampleStats.benchmarksharpe;
numbers(3,1)=SampleStats.mean_bl;
numbers(3,2)=SampleStats.mean_benchmarkport;
numbers(4,1)=SampleStats.max_bl;
numbers(4,2)=SampleStats.benchmarkmax;
numbers(5,1)=SampleStats.min_bl;
numbers(5,2)=SampleStats.benchmarkmin;

%Bl compared to Backwards Max Sharpe
numbers2 = zeros(5,2);

% Fill out the table with the data from SampleStats
numbers2(1,1)=SampleStats.cumret_bl(end,:);
numbers2(1,2)=SampleStats.cumret_ms(end,1);
numbers2(2,1)=SampleStats.sharpe_bl;
numbers2(2,2)=SampleStats.sharpe_ms;
numbers2(3,1)=SampleStats.mean_bl;
numbers2(3,2)=SampleStats.mean_ms;
numbers2(4,1)=SampleStats.max_bl;
numbers2(4,2)=SampleStats.max_ms;
numbers2(5,1)=SampleStats.min_bl;
numbers2(5,2)=SampleStats.min_ms;


%Compare additional risk metrics 
numbers3 = zeros(6,3);

% Fill out the table with the data from SampleStats
numbers3(1,1)=SampleStats.vol_bl;
numbers3(1,2)=SampleStats.vol_benchmarkport;
numbers3(1,3)=SampleStats.vol_ms;
numbers3(2,1)=SampleStats.skew_bl;
numbers3(2,2)=SampleStats.skew_benchmarkport;
numbers3(2,3)=SampleStats.skew_ms;
numbers3(3,1)=SampleStats.kurt_bl;
numbers3(3,2)=SampleStats.kurt_benchmarkport;
numbers3(3,3)=SampleStats.kurt_ms;
numbers3(4,1)=SampleStats.maxdraw_bl;
numbers3(4,2)=SampleStats.benchmarkmax_draw;
numbers3(4,3)=SampleStats.maxdraw_ms;
numbers3(5,1)=SampleStats.histVaR_bl;
numbers3(5,2)=SampleStats.histVaR_benchmarkport;
numbers3(5,3)=SampleStats.histVaR_ms;
numbers3(6,1)=SampleStats.gaussianVaR_bl;
numbers3(6,2)=SampleStats.gaussianVaR_benchmarkport;
numbers3(6,3)=SampleStats.gaussianVaR_ms;




%close all; % close figures

%Refer to SampleStats Structure for performance metrics. 

   %% Portfolio Construction   
    
    clearvars -except numbers numbers2 numbers3  Reflection  SampleStats startingCash student_name  hwgt hret hrsk swgt srsk sret cret crsk mret mrsk p historical holdingperiod macroCommentary marketCommentary mrkthistorical mrktholding PresName reportNumber symbols


%% Housekeeping
import mlreportgen.ppt.* % import Matlab package to create ppt presentation

%% Build Presentation

slides = Presentation(PresName); % create ppt presentation


%% Title Slide
slide0 = add(slides,'Title Slide'); % 'Title Slide' tells the type of slide

contents = find(slides,'Title'); % find Title object and assign to placeholder variable
replace(contents(1),reportNumber); % replace contents with text
contents(1).FontColor = 'blue';
replace(slide0,'Subtitle','Connor Haines');

%% Slide 1
slide1 = add(slides,'Title and Content');
replace(slide1,'Title','Portfolio Description')
replace(slide1,'Content','This is a long-only, buy-and-hold portfolio made up of US Mid Cap stocks. Portfolio optimized Sharpe Ratio by using the output from a Black-Litterman model with views formed through a combination of ARMA GARCH forecasting and discretionary perspectives. Portfolio was formed on 1/12/17 and has held all positions constant since then.');

%% Slide 2
slide2 = add(slides,'Title and Content');
replace(slide2,'Title','Allocation Strategy')
replace(slide2,'Content',{'Step 1: Use discretion to choose 20 possible US Mid Cap Stocks. All stocks chosen make up a material portion of IJH, the benchmark portfolio.','Step 2: Forecast returns using ARMA-GARCH forecasting.','Step 3: Use Black-Litterman model with views formed as a mixture of forecasted returns and discretionary views to optimize max sharpe ratio portfolio.', 'Step 4: Evaluate portfolio performance relative to the benchmark and a historical max sharpe ratio portfolio, paying particular attention to measures of risk.'})


%% Title Slide
slide3 = add(slides,'Title Slide');
replace(slide3,'Title','Market Commentary')

%% Slide 4
slide4 = add(slides,'Title and Content');
replace(slide4,'Title','Macro Commentary')
replace(slide4,'Content',macroCommentary);

%% Slide 5
slide5 = add(slides,'Title and Content');
replace(slide5,'Title','Market Commentary')
replace(slide5,'Content',marketCommentary);

%% Title Slide
slide6 = add(slides,'Title Slide');
replace(slide6,'Title','Portfolio-level Analysis')

%% Picture slide 2
close all; % To prevent picture from staying open

% Create slide and insert picture
Pie_alloc = Picture('Pie_alloc_bl.png');
picSlide2 = add(slides,'Title and Picture');
replace(picSlide2,'Title','Asset Allocation');
replace(picSlide2,'Picture',Pie_alloc);

%% Picture slide 3
% Create slide and insert picture
Hist_port = Picture('Hist_port.png');
picSlide3 = add(slides,'Title and Picture');
replace(picSlide3,'Title','Histogram of Returns');
replace(picSlide3,'Picture',Hist_port);

%% Picture slide 4
% Create slide and insert picture
Daily_ret = Picture('Daily_ret.png');
picSlide4 = add(slides,'Title and Picture');
replace(picSlide4,'Title','Daily Returns Compared to Benchmark');
replace(picSlide4,'Picture',Daily_ret);

%% Picture slide 5
% Create slide and insert picture
Cum_ret = Picture('Cum_ret.png');
picSlide5 = add(slides,'Title and Content');
replace(picSlide5,'Title','Cumulative Returns Compared to Benchmark');
replace(picSlide5,'Content',Cum_ret);

%% Table slide 1
% Create table
headers1 = {'Portfolio'; 'Benchmark'}';
headers2 = {'Since Inception';'Since Inception'}';
rownames = {'';'';'Cumulative Return';'Sharpe Ratio'; ...
            ;'Mean Return'; ...
            ;'Max Port Value';'Min Port Value';};
data1 = [headers1;headers2;num2cell(numbers(:,1:2))];
data = [rownames, data1];
myTable1 = Table(data);

% Create slide and insert table
tableSlide = add(slides,'Title and Content');
replace(tableSlide,'Title','Performance Statistics');
replace(tableSlide,'Content',myTable1); % remember to wrap the table in Table()


%% Title Slide
slide7 = add(slides,'Title Slide');
replace(slide7,'Title','Black Litterman Portfolio compared to Backwards Max Sharpe Portfolio')

%% Picture slide 6
close all; % To prevent picture from staying open

% Create slide and insert picture
Pie_alloc_ms = Picture('Pie_alloc_ms.png');
picSlide6 = add(slides,'Title and Picture');
replace(picSlide6,'Title','Asset Allocation - Black Litterman vs. Backwards Max Sharpe');
replace(picSlide6,'Picture',Pie_alloc_ms);


%% Picture slide 7
% Create slide and insert picture
Daily_ret_ms = Picture('Daily_ret_ms.png');
picSlide7 = add(slides,'Title and Picture');
replace(picSlide7,'Title','Daily Returns Compared to Backwards Max Sharpe');
replace(picSlide7,'Picture',Daily_ret_ms);

%% Picture slide 8
% Create slide and insert picture
Cum_ret_ms = Picture('Cum_ret_ms.png');
picSlide8 = add(slides,'Title and Content');
replace(picSlide8,'Title','Cumulative Returns Compared to Backwards Max Sharpe');
replace(picSlide8,'Content',Cum_ret_ms);


%% Table slide 2
% Create table
headers1 = {'Portfolio'; 'Backwards Max Sharpe'}';
headers2 = {'Since Inception';'Since Inception'}';
rownames = {'';'';'Cumulative Return';'Sharpe Ratio'; ...
            ;'Mean Return'; ...
            ;'Max Port Value';'Min Port Value';};
data1 = [headers1;headers2;num2cell(numbers2(:,1:2))];
data = [rownames, data1];
myTable2 = Table(data);

% Create slide and insert table
tableSlide = add(slides,'Title and Content');
replace(tableSlide,'Title','Performance Statistics');
replace(tableSlide,'Content',myTable2); % remember to wrap the table in Table()


%% Title Slide
slide8 = add(slides,'Title Slide');
replace(slide8,'Title','Risk Comparisson Across Portfolios')


%% Picture slide 9
% Create slide and insert picture
Hist_var = Picture('Hist_var.png');
picSlide9 = add(slides,'Title and Content');
replace(picSlide9,'Title','Historical Var - (top to bottom) Black Litterman - Benchmark - Historial Max Sharpe');
replace(picSlide9,'Content',Hist_var);

%% Table slide 3
% Create table
headers1 = {'Portfolio'; 'Benchmark'; 'Backwards Max Sharpe'}';
headers2 = {'Since Inception';'Since Inception';'Since Inception'}';
rownames = {'';'';'Standard Deviation';'Skew'; ...
            ;'Kurtosis'; ...
            ;'Max Drawdown';'Historical VaR';'Gaussian VaR';};
data1 = [headers1;headers2;num2cell(numbers3(:,1:3))];
data = [rownames, data1];
myTable3 = Table(data);

% Create slide and insert table
tableSlide = add(slides,'Title and Content');
replace(tableSlide,'Title','Risk Statistics');
replace(tableSlide,'Content',myTable3); % remember to wrap the table in Table()


%% Slide 9
slide9 = add(slides,'Title and Content');
q = Paragraph();
Reflection = Text('Of the three portfolios, my portfolio performed the best, outperforming the benchmark and historical max sharpe ratio portfolio in terms of cumulative returns and sharpe ratio. Overtime, my portfolio pulled away from the benchmark portfolio and ended up with close to 6% higher returns. However, throughout the investment horizon, my portfolio tracked closely with the historical max sharpe ratio portfolio. I think this makes sense. Both the black litterman model and the historical sharpe ratio portfolio use a variation of finding a historical efficient frontier. BL does this through starting with a market equilibrium distribution of assets, and then using a historical CAPM model (which has an efficient frontier implicit within it) combined with the user inputed views. This returns a BL output that is then inputted into a sharpe ratio optimizer. The real advantage to using a Black Litterman model is that it results in a less concentrated asset distribution than the historical max sharpe ratio portfolio. The asset distributions on slide 14 show this clearly. The risk analysis on the previous slide shows that so far, the risk characteristics of the two portfolios have been similar, however over a longer time horizon and/or if markets become more volatile, the greater diversification in the Black Literrman model could prove useful. Additionally, views inputted into the Black Litterman framework have resulted in a portfolio that has outperformed the historical max sharpe ratio portfolio by ~1.7% since the beginning of the investment horizon.');
Reflection.Style = {FontSize('18pt')};
append(q,Reflection);
replace(slide9,'Content', q);
replace(slide9,'Title','Closing Thoughts');

%% Generate and open the presentation
close(slides);

if ispc % if this machine is Windows
    winopen(PresName);
end


clearvars -except SampleStats startingCash hwgt hret hrsk swgt srsk sret cret crsk mret mrsk p pret prsk pwgt Ret historical holdingperiod  mrkthistorical mrktholding PresName reportNumber symbols








    
   