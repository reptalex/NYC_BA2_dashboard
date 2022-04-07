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
NYC[,raw_cases:=new_confirmed] ### raw data
NYC[,new_confirmed:=outlier_detection(new_confirmed)]
NYC[,rm:=round(frollmean(new_confirmed,7,align='right'))]

### Missing data on 4/3 and a next-day data dump on 4/4 are outliers
NYC[date %in% as.Date(c('2022-04-03','2022-04-04')),new_confirmed:=NA]
NYC[date %in% as.Date('2022-04-03'),rm:=NA]

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

ggplot(xx[outbreak_time > -10],aes(outbreak_time,growth_rate))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=p2.5_growth_rate,ymax=p97.5_growth_rate,fill=region),alpha=0.1)+
  geom_line(aes(color=region),lwd=2)+
  scale_color_manual(values=c('black','red'))+
  scale_fill_manual(values=c('black','red'))+
  ggtitle('UK & NYC BA.2 Outbreaks, case growth rates')+
  scale_y_continuous('Case exponential growth rate, r(t)')+
  scale_x_continuous('Outbreak Time (days since r(t)>0)')+
  theme(legend.position=c(0.15,0.9))
ggsave('figures/Outbreak_daily_growth_rate_comparison.png',height=8,width=10)

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
#   theme(legend.position=c(0.15,0.9))+
#   geom_line(data=nyc_out_of_sample,lty=2,col='red',lwd=2)
# grs
# ggsave('figures/Outbreak_7dma_growth_rate_comparison.png',height=8,width=10)
# 
# 

# forecasting NYC ---------------------------------------------------------
nyc <- xx[region=='New York City']
uk <- xx[region=='United Kingdom']

nyc <- forecast(nyc,uk) %>% glue_ma

medium_alert_threshold=(200/1e5)*8.4e6/7


ll <- function(sd,
               mn=log(nyc[.N,rm]),
               up=log(nyc[.N,p97.5_rm]),
               lw=log(nyc[.N,p2.5_rm])){
  ((up-qnorm(.975,mean=mn,sd=sd))^2+(lw-qnorm(.025,mean=mn,sd=sd))^2) %>% return
}

sdlog <- optimize(f=ll,interval = c(.1,1e3))$minimum

prob=signif(1-plnorm(medium_alert_threshold,meanlog=log(nyc[.N,rm]),sdlog=sdlog),2)*100


ggplot(nyc[forecast==TRUE],aes(date,exp(mean_position)))+
  geom_line(data=nyc[forecast==FALSE],lwd=2,col='red')+
  geom_line(lwd=2,col='red',lty=2)


fc <- ggplot(nyc,aes(date,rm))+
  geom_line(data=nyc[forecast==FALSE & !is.na(rm)],lwd=2,col='red')+
  geom_line(data=nyc[forecast==TRUE],lwd=2,lty=2,col='red')+
  # geom_bar(data=nyc,aes(y=new_confirmed),stat='identity',col='red',fill=rgb(1,0,0,0.3))+
  geom_ribbon(data=nyc[forecast==TRUE],aes(ymin=p2.5_rm,ymax=p97.5_rm),fill='red',alpha=0.2)+
  geom_hline(yintercept = medium_alert_threshold,lwd=1)+
  annotate(geom='text',x=as.Date('2022-03-18'),y=medium_alert_threshold+90,
           label="Medium Alert Level Threshold",size=6)+
  scale_y_continuous('New Cases')+
  annotate(geom='segment',x=as.Date('2022-04-03'),xend=max(nyc_cases$date),y=2620,yend=medium_alert_threshold)+
  annotate(geom='text',x=as.Date('2022-04-03'),y=2730,
           label=paste(prob,'% chance of\n Medium Alert by ',max(nyc_cases$date),sep=''),
           fill='white',color='black')+
  geom_point(data=ny_tot_manual,cex=4,color='darkred')+
  geom_line(data=ny_tot_manual,color='darkred')+
  ggtitle('NYC Outbreak Forecast: If NYC follows UK trajectory')

ggarrange(grs,fc,ncol=2,labels=c("A",'B'))
ggsave('figures/BA2_NYC_forecast.png',height=6,width=15)




nyc_cases[p97.5>medium_alert_threshold,min(date)]
