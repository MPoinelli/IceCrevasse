% Mattia Poinelli
% JPL, April 2019
% Summary of LEFM results and elaboration of features' characteristics table
% Taking into account cycloid 1-2-4-6
% SI units, unless differently specified!!!

format short

MAX_PRO_RATE = [max(CRACK(1).PROPAGATION_RATES),max(CRACK(2).PROPAGATION_RATES),...
    max(CRACK(4).PROPAGATION_RATES),max(CRACK(6).PROPAGATION_RATES)]';

MAX_OPE_RATE = [max(CRACK(1).OPENING_RATE),max(CRACK(2).OPENING_RATE),...
    max(CRACK(4).OPENING_RATE),max(CRACK(6).OPENING_RATE)]';

OPENING_WIDT = [max(CRACK(1).WIDTH),max(CRACK(2).WIDTH),...
    max(CRACK(4).WIDTH),max(CRACK(6).WIDTH)]';

TIME_TO_DEVE = [CRACK(1).TIME(end-1),CRACK(2).TIME(end-1),CRACK(4).TIME(end-1),CRACK(6).TIME(end-1)]';

STANDBY_PERI = [sum(CRACK(1).STANDBY),sum(CRACK(2).STANDBY),sum(CRACK(4).STANDBY),sum(CRACK(6).STANDBY)]';

LENGTH       = [CRACK(1).LENGTH,CRACK(2).LENGTH,CRACK(4).LENGTH,CRACK(6).LENGTH]';

GROWTH_RATE  = LENGTH./TIME_TO_DEVE;

% conversion
MAX_OPE_RATE_mm_yr = MAX_OPE_RATE.*(3600*24*365*1000);
TIME_TO_DEVE_days  = TIME_TO_DEVE./(3600*24);
STANDBY_PERI_days  = STANDBY_PERI./(3600*24);
GROWTH_RATE_m_day  = GROWTH_RATE.*(36*24);

T = table(MAX_PRO_RATE,MAX_OPE_RATE_mm_yr,TIME_TO_DEVE_days,OPENING_WIDT,STANDBY_PERI_days,GROWTH_RATE_m_day,'RowNames',{'1','2','Delphi','Cilicia'})