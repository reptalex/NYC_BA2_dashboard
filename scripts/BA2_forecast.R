library(tidyverse)
library(COVID19)
library(padr)
library(lubridate)
library(progress)
library(RcppRoll)
library(parallel)
library(data.table)
library(ggpubr)
library(usmap)
library(tsoutliers)
library(mgcv)
library(gganimate)
source("scripts/utils.R")


# Settings & params -------------------------------------------------------

theme_set(theme_bw(base_size=12))
min_date <- as.Date('2020-10-01')
filt=TRUE
forecast_after <- as.Date('2022-03-31')


# Countries ---------------------------------------------------------------

UK <- covid19('United Kingdom') %>% as.data.table
UK[,country:=administrative_area_level_1]
UK[,region:=country]
UK <- UK[country %in% countries & date>=min_date]
UK[,raw_cases:=dfs(confirmed),by=country]
UK[,new_confirmed:=outlier_detection(raw_cases)]
UK[,rm:=round(frollmean(new_confirmed,7,align='right'))]

UK[,new_deaths:=dfs(deaths),by=country]
UK[,deaths_pc:=deaths/population]

UK <- covid19_nbss(UK,filtering = filt) %>% as.data.table

# US holidays -------------------------------------------------------------

thanksgiving <- c(seq(as.Date('2020-11-26'),as.Date('2020-11-27'),by='day'),
                  seq(as.Date('2021-11-25'),as.Date('2021-11-29'),by='day'))
christmas <- as.Date(c('2020-12-25','2020-12-26'))
ny<- seq(as.Date('2021-01-01'),as.Date('2021-01-02'),by='day')
lbr <- as.Date(c('2021-09-06','2021-09-07'))
ipd <- as.Date(c('2021-10-11','2021-10-12'))
hlwn <- as.Date(c('2021-10-30','2021-10-31','2021-11-01'))
us_holidays <- c(thanksgiving,christmas,ny,lbr,ipd,hlwn)

# NYC ----------------------------------------------------
US_counties <- covid19("United States",level=3) %>% as.data.table
US_counties[,country:=administrative_area_level_1]
US_counties[,state:=administrative_area_level_2]
US_counties[,county:=administrative_area_level_3]
US_counties <- US_counties[!state %in% c('Guam','Northern Mariana Islands','Virgin Islands','American Samoa')]
US_counties[,new_confirmed:=dfs(confirmed),by=c('state','county')]

US_counties[,new_deaths:=dfs(deaths),by=c('state','county')]
US_counties[,lbl:=paste(state,county,sep=',')]

NYC <- US_counties[county=='New York City' & date>=min_date]

### there's an important reporting error in the NYC data.
### On Sunday, April 3rd, NYC broke from a previous pattern of reporting on Sundays & reported 0 cases
### The following day, April 4th, NYC had a Sunday-like spike of reporting. 
### This outlier of changing day-of-week reporting comes at a critical point of the NYC outbreak
### when decelerations were first evident and a peak could occur in the following week or two
### One approach would be to forecast 7d moving averages...
### but such moving average forecasts lag behind contemporary incidence, underestimating incidence
### at the peak and overestimating incidence after the peak.
### To keep forecasts up-to-date and accurately incorporating new daily information with minimal lag,
### I will manually split the Monday April 4th cases across Sunday April 3 and Monday April 4,
### preserving the Sun/Mon case ratio estimated from the preceding 4 weeks.
NYC[,raw_cases:=new_confirmed] ### raw data - for later use if needed.

Sunday_cases <- NYC[date %in% (as.Date('2022-04-03')-7*(1:4)),new_confirmed]
Monday_cases <- NYC[date %in% (as.Date('2022-04-04')-7*(1:4)),new_confirmed]
Anomalous_date_cases <- NYC[date==as.Date('2022-04-04'),new_confirmed]
Sun_Mon_ratio <- exp(mean(log(Sunday_cases/Monday_cases)))

NYC[date==as.Date('2022-04-03')]$new_confirmed <- round(Anomalous_date_cases*Sun_Mon_ratio/(Sun_Mon_ratio+1))
NYC[date==as.Date('2022-04-04')]$new_confirmed <- round(Anomalous_date_cases/(Sun_Mon_ratio+1))

NYC[,new_confirmed:=outlier_detection(new_confirmed)]

NYC[,rm:=round(frollmean(new_confirmed,7,align='right'))]

### Missing data on 4/3 and a next-day data dump on 4/4 are outliers
# NYC[date %in% as.Date(c('2022-04-03','2022-04-04')),new_confirmed:=NA]
# NYC[date %in% as.Date('2022-04-03'),rm:=NA]

NYC <- covid19_nbss(NYC,filtering = filt) %>% as.data.table


# Combining datasets ----------------------------------------------------------------

uk <- UK[date>as.Date('2022-02-01')]
uk[,region:='United Kingdom']
nyc <- NYC[date>as.Date('2022-02-01')]
nyc[,region:='New York City']

nyc_tot <- nyc
# nyc <- nyc[date<=forecast_after]

xx <- rbind(uk[,c('date','region','growth_rate','p2.5_growth_rate','p97.5_growth_rate',
                  'new_confirmed','rm','raw_cases','mean_position','p2.5_position','p97.5_position')],
            nyc[,c('date','region','growth_rate','p2.5_growth_rate','p97.5_growth_rate',
                   'new_confirmed','rm','raw_cases','mean_position','p2.5_position','p97.5_position')])
xx[,region:=factor(region,levels=c('United Kingdom','New York City'))]



# Outbreak comparison -----------------------------------------------------

xx[,outbreak_start:=date[min(which(growth_rate>0 & date>as.Date('2022-02-14')))],by=region]
xx[,outbreak_time:=as.numeric(date-outbreak_start)]

gr <- ggplot(xx[outbreak_time > -10],aes(outbreak_time,growth_rate))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=p2.5_growth_rate,ymax=p97.5_growth_rate,fill=region),alpha=0.1)+
  geom_line(aes(color=region),lwd=2)+
  scale_color_manual(values=c('black','red'))+
  scale_fill_manual(values=c('black','red'))+
  ggtitle('UK & NYC BA.2 Outbreaks, case growth rates')+
  scale_y_continuous('Case exponential growth rate, r(t)')+
  scale_x_continuous('Outbreak Time (days since r(t)>0)')+
  theme(legend.position=c(0.15,0.9))
gr
ggsave('figures/Outbreak_daily_growth_rate_comparison.png',height=8,width=10)



# Outbreak comparison GIF -------------------------------------------------
p=ggplot(xx[outbreak_time> -10],
         aes(outbreak_time,growth_rate,group=region,frame=date,fill=region))+
  transition_reveal(date)+
  geom_line(lwd=2,aes(color=region))+
  geom_ribbon(aes(ymin=p2.5_growth_rate,ymax=p97.5_growth_rate,fill=region),alpha=0.2)+
  enter_fade()+
  scale_x_continuous('Outbreak Time',limits=c(-10,35))+
  exit_fade()+
  scale_y_continuous('Growth Rate',limits=c(-0.15,0.15))+
  ggtitle('BA.2 Outbreak Tracker | Date: {frame_along}')+
  scale_color_manual(values=c('black','red'))+
  scale_fill_manual(values=c('black','red'))+
  theme(legend.position=c(0.8,0.8))+
  geom_hline(yintercept=0,color='darkgrey')

animate(p,nframes = length(unique(xx$date)),fps = 8) %>%
  anim_save(filename='figures/UK_NYC_BA_2_outbreaks.gif',animation=.)


# 7d moving average growth rate estimation --------------------------------
### Due to changing day-of-week reporting patterns
### we can't use the default day-of-week adjusted growth rate estimator.
### Instead, we focus on forecasting the growth rate of the 7-day moving average.

# xx[is.na(new_confirmed),new_confirmed:=0]
# xx[region=="New York City" & date==as.Date('2022-04-03'),rm:=NA] ### an unusual Sunday without reported cases
# xx[,grm:=nbs(rm,filtering=TRUE,seasonal=FALSE),by=region]
# xx[,grm_2.5:=nbs(rm,filtering=TRUE,seasonal=FALSE,name = 'p2.5_growth_rate'),by=region]
# xx[,grm_97.5:=nbs(rm,filtering=TRUE,seasonal=FALSE,name = 'p97.5_growth_rate'),by=region]
# 
# xx[,outbreak_start_date:=date[min(which(grm>0 & date>as.Date('2022-02-20')))],by=region]
# xx[,outbreak_time:=as.numeric((date-outbreak_start_date))]

# nyc_out_of_sample[,growth]

# Outbreak r(t) comparison ------------------------------------------------------------
# 
# 
# grs <- ggplot(xx[outbreak_time > -10],aes(outbreak_time,grm))+
#   geom_hline(yintercept = 0)+
#   geom_ribbon(aes(ymin=grm_2.5,ymax=grm_97.5,fill=region),alpha=0.1)+
#   geom_line(aes(color=region),lwd=2)+
#   scale_color_manual(values=c('black','red'))+
#   scale_fill_manual(values=c('black','red'))+
#   ggtitle('UK & NYC BA.2 Outbreaks, 7d ma growth rates')+
#   scale_y_continuous('Case exponential growth rate, r(t)')+
#   scale_x_continuous('Outbreak Time (days since r(t)>0)')+
#   theme(legend.position=c(0.15,0.9))
# grs
# ggsave('figures/Outbreak_7dma_growth_rate_comparison.png',height=8,width=10)
# 
# 

# forecasting NYC ---------------------------------------------------------
nyc <- xx[region=='New York City']
uk <- xx[region=='United Kingdom']

nyc <- forecast(nyc,uk,smooth_method='GLM') %>% glue_ma

medium_alert_threshold=(200/1e5)*8.4e6/7


ll <- function(sd,
               mn=log(nyc[.N,rm]),
               up=log(nyc[.N,p97.5_rm]),
               lw=log(nyc[.N,p2.5_rm])){
  ((up-qnorm(.975,mean=mn,sd=sd))^2+(lw-qnorm(.025,mean=mn,sd=sd))^2) %>% return
}

sdlog <- optimize(f=ll,interval = c(.1,1e3))$minimum
prob=signif(1-plnorm(medium_alert_threshold,meanlog=log(nyc[.N,rm]),sdlog=sdlog),2)*100


# Forecast figure ---------------------------------------------------------

### This figure requires some manual curation to format dates & annotate the plot. 
nyc[,Date:=date] ### changes X-axis to "Date" instead of "date" w/o upsetting the x scale
mxd=max(nyc$date) ### print to upate max_date_char, and use this to label confidence intervals
mxd
max_date_char <- "April 19, 2022"
cd <- min(nyc[forecast==TRUE,date])
x1 <- nyc[forecast==FALSE,rm[.N]]
mn_ci <- nyc[.N-1,p2.5_rm]
fc <- ggplot(nyc[forecast==TRUE],aes(Date,rm))+
  geom_line(data=nyc[forecast==FALSE & !is.na(rm) & outbreak_time>-10],lwd=2,col='red')+
  geom_line(lwd=2,lty=2,col='red')+
  geom_bar(data=nyc[forecast==FALSE & outbreak_time > -10],aes(y=new_confirmed),col='red',fill='red',alpha=0.1,stat='identity')+
  # geom_bar(data=nyc,aes(y=new_confirmed),stat='identity',col='red',fill=rgb(1,0,0,0.3))+
  geom_ribbon(data=nyc[forecast==TRUE],aes(ymin=p2.5_rm,ymax=p97.5_rm),fill='red',alpha=0.2)+
  geom_hline(yintercept = medium_alert_threshold,lwd=1)+
  annotate(geom='text',x=as.Date('2022-03-18'),y=medium_alert_threshold+90,
           label="Medium Alert Level Threshold",size=6)+
  scale_y_continuous('Newly Reported Cases',limits=c(0,3900))+
  # scale_x_continuous('Date')+
  annotate(geom='segment',x=as.Date('2022-04-10'),xend=max(nyc$date),y=3120,yend=medium_alert_threshold)+
  annotate(geom='text',x=as.Date('2022-04-10'),y=3290,
           label=paste(prob,'% chance of reaching \n "medium alert" by ',max_date_char,sep=''),
           fill='white',color='black')+
  annotate(geom='segment',x=cd,xend=cd+4,y=x1-20,yend=1200,color='red',lwd=2)+
  annotate(geom='text',x=cd+5,y=1100,
           label='7-day moving average',
           fill='red',color='red',bg='red')+
  annotate(geom='segment',x=mxd-1,xend=mxd-2,y=mn_ci,yend=800,color='red',lwd=2,alpha=0.2)+
  annotate(geom='text',x=mxd-2,y=750,
           label='95% Confidence Intervals',
           fill='red',color='red',bg='red')+
  theme_classic(base_size=15)+
  # geom_point(data=ny_tot_manual,cex=4,color='darkred')+
  # geom_line(data=ny_tot_manual,color='darkred')+
  ggtitle('Forecast of New York City COVID-19 Cases over next 12 days')

fc
ggsave('figures/BA2_NYC_forecast.png',height=8,width=11)


gr+
  geom_line(data=nyc[forecast==TRUE],lwd=2,col='red',lty=2)
ggsave('figures/Outbreak_daily_growth_rate_comparison_with_forecast.png',height=6,width=8)

save(list=ls(),file='data/BA2_forecast_workspace.Rds')
