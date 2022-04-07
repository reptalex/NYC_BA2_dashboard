### Utils from https://github.com/reptalex/COVID_ERL, with modified outlier detection


### These tools were used in the following papers:
## Analysis of COVID outbreaks on a timescale of burden:
## https://www.medrxiv.org/content/10.1101/2021.05.03.21256542v2

## Characaterizing features of outbreak duration for novel SARS-CoV-2 variants of concern
## https://www.medrxiv.org/content/10.1101/2022.01.14.22269288v1


library(data.table)
library(deSolve)
library(magrittr)
library(zoo)
library(KFAS)
library(parallel)
library(tidyverse)
library(tsoutliers)
library(progress)
library(lubridate)
library(EpiEstim)
library(mgcv)

dfs <- function(x){
  x <- c(x[1],diff(x))
  x[x<0] <- NA
  return(x)
}

# Requires to get sockets parallization to work on os x. 
library(rstudioapi)
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  if(versionInfo()$version < "1.3.1056"){
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  }  
}

nbs <- function(x,name='growth_rate',seasonal=TRUE,remove_initial_zeros=TRUE,...){
  y <- rep(NA,length(x))
  if (remove_initial_zeros){
    first_nz <- min(which(x>0))
  } else {
    first_nz <- 1
  }
  
  if (name=='all'){
    mdl <- nbss(x[first_nz:length(x)],...)
    
    dd <- matrix(NA,nrow=first_nz-1,ncol=ncol(mdl))
    colnames(dd) <- colnames(mdl)
    y <- rbind(dd,mdl)
  } else {
    y[first_nz:length(x)] <- nbss(x[first_nz:length(x)],...)  %>% getElement(name)
  }
  return(y)
}

nbss <- function(x,filtering=FALSE,seasonal=TRUE,dispersion=NULL){
  if (seasonal){
    nb_model <- function(x, pars){
      model_nb <- SSModel(x ~ SSMtrend(2, Q=list(0, NA),
                                       P1=diag(c(10, 1)),
                                       a1=c(0, 0),
                                       state_names=c("level", "trend"))+
                            SSMseasonal(7),
                          u=rep(exp(pars[1]), length(x)), distribution="negative binomial")
      fit <- fitSSM(model_nb, c(0), method="L-BFGS-B", control=list(maxit=200))
      return(fit)
    }
  } else {
    nb_model <- function(x, pars){
      model_nb <- SSModel(x ~ SSMtrend(2, Q=list(0, NA),
                                       P1=diag(c(10, 1)),
                                       a1=c(0, 0),
                                       state_names=c("level", "trend")),
                          u=rep(exp(pars[1]), length(x)), distribution="negative binomial")
      fit <- fitSSM(model_nb, c(0), method="L-BFGS-B", control=list(maxit=200))
      return(fit)
    }
  }
  
  if (is.null(dispersion)){
    logLik_nb <- function(x, pars){
      fit <- nb_model(x, pars)
      ll <- logLik(fit$model, marginal = TRUE)
      return(-ll)
    }
    
    res <- tryCatch(optim(c(-1), function(y) logLik_nb(x, y), method="Brent", lower=-2, upper=2) , error=function(e) NULL)
    if (is.null(res)) return(NULL)
    fit <- nb_model(x, res$par)
  } else {
    fit <- tryCatch(nb_model(x,dispersion),error=function(e) NULL)
    res <- NULL
    res$par[1] <- dispersion
  }
  if (filtering==TRUE){
    sm_signal <- KFS(fit$model, filtering="signal",smoothing='none')
    sm_state <- KFS(fit$model, filtering="state",smoothing='none')
    
    out <- data.frame(p2.5_position = c(qnorm(0.025, sm_state$a[-1,'level'], (sqrt(sm_state$P[1,1,-1])))), 
                      p97.5_position = c(qnorm(0.975, sm_state$a[-1,'level'],(sqrt(sm_state$P[1,1,-1])))), 
                      mean_position = (c(sm_state$a[-1,'level'])),
                      p2.5_growth_rate = c(qnorm(0.025, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                      p97.5_growth_rate = c(qnorm(0.975, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                      p25_growth_rate = c(qnorm(0.25, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                      p75_growth_rate = c(qnorm(0.75, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                      growth_rate = (c(sm_state$a[-1,'trend'])),
                      percentile_0_growth_rate =c(pnorm(0, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                      growth_rate = (c(sm_state$a[-1,'trend'])),
                      dispersion=res$par[1], 
                      z_score_growth_rate = c(sm_state$a[-1,2]/sqrt(sm_state$P[2,2,-1])))
  } else {
    sm_signal <- KFS(fit$model, smoothing="signal")
    sm_state <- KFS(fit$model, smoothing="state")
    out <- data.frame(p2.5_signal = exp(qnorm(0.025, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                      p97.5_signal = exp(qnorm(0.975, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                      mean_signal = exp(sm_signal$thetahat), 
                      p2.5_position = c(qnorm(0.025, sm_state$alphahat[,1], (sqrt(sm_state$V[1,1,])))), 
                      p97.5_position = c(qnorm(0.975, sm_state$alphahat[,1],(sqrt(sm_state$V[1,1,])))), 
                      mean_position = (c(sm_state$alphahat[,1])),
                      p2.5_growth_rate = c(qnorm(0.025, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                      p97.5_growth_rate = c(qnorm(0.975, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                      p25_growth_rate = c(qnorm(0.25, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                      p75_growth_rate = c(qnorm(0.75, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                      growth_rate = (c(sm_state$alphahat[,2])),
                      percentile_0_growth_rate =c(pnorm(0, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                      growth_rate = (c(sm_state$alphahat[,2])),
                      dispersion=res$par[1], 
                      z_score_growth_rate = c(sm_state$alphahat[,2]/sqrt(sm_state$V[2,2,])))
  }
  return(out)
}


fit_covid_ssm <- function(d, series="new_confirmed", precomputed_dispersions=NULL,
                          return_fit=FALSE,maxiter=200,filtering=FALSE){
  dat <- d
  
  # Helper functions
  nb_model <- function(dat, pars){
    model_nb <- SSModel(dat[,series] ~ SSMtrend(2, Q=list(0, NA),
                                                     P1=diag(c(10, 1)),
                                                     a1=c(0, 0),
                                                     state_names=c("level", "trend"))+
                          SSMseasonal(7),
                        u=rep(exp(pars[1]), nrow(dat)), distribution="negative binomial")
    fit <- fitSSM(model_nb, c(0), method="L-BFGS-B", control=list(maxit=200))
    return(fit)
  }
  logLik_nb <- function(dat, pars){
    fit <- nb_model(dat, pars)
    ll <- logLik(fit$model, marginal = TRUE)
    return(-ll)
  }
  
  # Remove preceding zeros
  dat <- dat %>% 
    arrange(date)
  dat$cs = cumsum(ifelse(is.na(dat[,series]), 0, dat[,series]))
  pass <- dat$cs > 0 # rows to keep in analysis 
  # dat <- dat %>% filter(cs > 0)
  
  
  if (nrow(dat[pass,]) < 10) return(cbind(d, error="not enough non-zero"))
  
  
  # Remove weekend zeros
  dat[[series]] <- ifelse((dat[[series]] == 0) & (format(dat$date, "%u") %in% c(6,7)), 
                          NA, dat[[series]])
  # Remove holiday zeros
  dat[[series]] <- ifelse((dat[[series]] == 0) & (dat$date %in% c(ymd("2020-11-26"), # Thanksgiving
                                                                  ymd("2020-12-25"), 
                                                                  ymd("2020-12-31"), 
                                                                  ymd("2021-01-01"))),
                          NA, dat[[series]])
  
  
  # outlier detection for early outbreak
  pass <- custom_processors(dat, pass)
  if (sum(dat[pass,series]!=0, na.rm=TRUE) < 10) return(cbind(d, error="not enough non-zero"))
  tryCatch({
    tmp <-getElement(dat,series)[pass]
    #tmp <- ifelse(is.na(tmp), 0, tmp)
    tmp <- na.approx(tmp)
    outlier_filtered_ts <- tmp %>% outlier_detection
    outlier_filtered_ts[is.na(getElement(dat, series)[pass])] <- NA
    tmp <- rep(NA, nrow(dat))
    tmp[pass] <- outlier_filtered_ts
    dat <- mutate(dat, series=tmp)    
  },  error = function(err){
    return(cbind(d, error="outlier detection errored"))
  })
  
  if (length(dat[pass,series][!is.na(dat[pass,series])]) < 10) return(cbind(d, error="too many NA"))
  
  
  if (!is.null(precomputed_dispersions)){
    res <- list()
    res$par <- dplyr::filter(precomputed_dispersions, id==unique(dat$id))$dispersion
    if (length(res$par)==0) return(cbind(d, error="could not find precomputed dispersions")) # no precomputed dispersion (previously was not able to fit likely)
  } else {
    res <- tryCatch(optim(c(-1), function(x) logLik_nb(dat[pass,], x), method="Brent", lower=-2, upper=3), 
                    error = function(e) NULL)
    if (is.null(res)) return(cbind(d, error="dispersion optimization failed"))
    if (res$convergence != 0) return(cbind(d, error="dispersion optimization failed"))
  }
  
  # now fit the model with the optimized dispersion parameters
  fit <- tryCatch(nb_model(dat[pass,], res$par),error=function(e) NA)
  if (is.na(fit)){
    return(NULL)
  } else {
    if(return_fit) return(fit)
    if (fit$optim.out$convergence != 0) return(cbind(d, error="model optimiztion (not dispersion) failed"))
    if (filtering==FALSE){
      sm_signal <- KFS(fit$model, smoothing="signal")
      sm_state <- KFS(fit$model, smoothing="state")
      
      if (series=="new_confirmed"){
        out <- data.frame(p2.5_signal = exp(qnorm(0.025, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                          p97.5_signal = exp(qnorm(0.975, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                          mean_signal = exp(sm_signal$thetahat), 
                          p2.5_position = c(qnorm(0.025, sm_state$alphahat[,1], (sqrt(sm_state$V[1,1,])))), 
                          p97.5_position = c(qnorm(0.975, sm_state$alphahat[,1],(sqrt(sm_state$V[1,1,])))), 
                          mean_position = (c(sm_state$alphahat[,1])),
                          p2.5_growth_rate = c(qnorm(0.025, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                          p97.5_growth_rate = c(qnorm(0.975, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                          p25_growth_rate = c(qnorm(0.25, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                          p75_growth_rate = c(qnorm(0.75, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                          growth_rate = (c(sm_state$alphahat[,2])),
                          percentile_0_growth_rate =c(pnorm(0, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                          dispersion=res$par[1], 
                          z_score_growth_rate = c(sm_state$alphahat[,2]/sqrt(sm_state$V[2,2,])))
      } else if (series == "new_deaths"){
        out <- data.frame(p2.5_signal_deaths = exp(qnorm(0.025, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                          p97.5_signal_deaths = exp(qnorm(0.975, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                          mean_signal_deaths = exp(sm_signal$thetahat), 
                          p2.5_position_deaths = c(qnorm(0.025, sm_state$alphahat[,1], (sqrt(sm_state$V[1,1,])))), 
                          p97.5_position_deaths = c(qnorm(0.975, sm_state$alphahat[,1],(sqrt(sm_state$V[1,1,])))), 
                          mean_position_deaths = (c(sm_state$alphahat[,1])),
                          p2.5_growth_rate_deaths = c(qnorm(0.025, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                          p97.5_growth_rate_deaths = c(qnorm(0.975, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                          p25_growth_rate_deaths = c(qnorm(0.25, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                          p75_growth_rate_deaths = c(qnorm(0.75, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                          growth_rate_deaths = (c(sm_state$alphahat[,2])),
                          percentile_0_growth_rate_deaths =c(pnorm(0, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                          dispersion_deaths=res$par[1], 
                          z_score_growth_rate_deaths = c(sm_state$alphahat[,2]/sqrt(sm_state$V[2,2,])))
      }
    } else {
      if (series=="new_confirmed"){
        sm_signal <- KFS(fit$model, filtering="signal",smoothing='none')
        sm_state <- KFS(fit$model, filtering="state",smoothing='none')
        
        out <- data.frame(p2.5_position = c(qnorm(0.025, sm_state$a[-1,'level'], (sqrt(sm_state$P[1,1,-1])))), 
                          p97.5_position = c(qnorm(0.975, sm_state$a[-1,'level'],(sqrt(sm_state$P[1,1,-1])))), 
                          mean_position = (c(sm_state$a[-1,'level'])),
                          p2.5_growth_rate = c(qnorm(0.025, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                          p97.5_growth_rate = c(qnorm(0.975, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                          p25_growth_rate = c(qnorm(0.25, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                          p75_growth_rate = c(qnorm(0.75, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                          growth_rate = (c(sm_state$a[-1,'trend'])),
                          percentile_0_growth_rate =c(pnorm(0, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                          growth_rate = (c(sm_state$a[-1,'trend'])),
                          dispersion=res$par[1], 
                          z_score_growth_rate = c(sm_state$a[-1,2]/sqrt(sm_state$P[2,2,-1])))
      } else if (series == "new_deaths"){
        
        out <- data.frame(p2.5_position_deaths = c(qnorm(0.025, sm_state$a[-1,'level'], (sqrt(sm_state$P[1,1,-1])))), 
                          p97.5_position_deaths = c(qnorm(0.975, sm_state$a[-1,'level'],(sqrt(sm_state$P[1,1,-1])))), 
                          mean_position_deaths = (c(sm_state$a[-1,'level'])),
                          p2.5_growth_rate_deaths = c(qnorm(0.025, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                          p97.5_growth_rate_deaths = c(qnorm(0.975, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                          p25_growth_rate_deaths = c(qnorm(0.25, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                          p75_growth_rate_deaths = c(qnorm(0.75, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                          growth_rate_deaths = (c(sm_state$a[-1,'trend'])),
                          percentile_0_growth_rate_deaths =c(pnorm(0, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                          growth_rate_deaths = (c(sm_state$a[-1,'trend'])),
                          dispersion_deaths=res$par[1], 
                          z_score_growth_rate_deaths = c(sm_state$a[-1,2]/sqrt(sm_state$P[2,2,-1])))
      }
    }
    
    if (any(grepl('administrative_area',colnames(dat)))){
      if (any(out$p97.5_signal > 1e8)) return(cbind(d, error="97.5_signal > 1e8"))
    }
    if (quantile(abs(out$z_score_growth_rate), probs=0.75) < 0.4) return(cbind(d, error="quantile(abs(out$z_score_growth_rate), probs=0.75) < 0.4"))
    
    tack_on <- matrix(NA, nrow=nrow(dat), ncol=ncol(out)) %>% as.data.frame()
    colnames(tack_on) <- colnames(out)
    tack_on[pass,] <- out
    return(cbind(dat, tack_on))
  }
}


covid19_nbss <- function(dat,series="new_confirmed", level='all',
                         mc.cores=1, precomputed_dispersions=NULL,
                         filtering=FALSE){
  if (level=="country"){
    tmp <- filter(dat, administrative_area_level==1)
  } else if (level=="state"){
    tmp <- filter(dat, administrative_area_level==2)
  } else if (level=="all") {
    tmp <- dat
  } else {
    stop("only level variables that are supported are all, country, and state")
  }
  tmp <- dat %>% 
    mutate(new_deaths = ifelse(new_deaths < 0, 0, new_deaths), 
           new_confirmed = ifelse(new_confirmed < 0, 0, new_confirmed)) %>% 
    as.data.frame() %>% 
    split(.$id)
  if (mc.cores == 1){
    fits <- list()
    pb <- progress_bar$new(total = length(tmp), format=" [:bar] :percent eta: :eta")
    for (i in 1:length(tmp)){
      pb$tick()
      fits[[i]] <- fit_covid_ssm(tmp[[i]], series, precomputed_dispersions,filtering=filtering)
      #if (is.null(fits[[i]])) stop("foo")
    }  
  } else {
    cl <- parallel::makeCluster(mc.cores)
    parallel::clusterEvalQ(cl, {
      library(tidyverse)
      library(lubridate)
      library(KFAS)
      library(tsoutliers)
      library(data.table)
    })
    parallel::clusterExport(cl,c("custom_processors", "outlier_detection", "fit_covid_ssm"))
    fits <- parLapply(cl, tmp,  function(x,series,precomputed_dispersions,filtering) fit_covid_ssm(x, series, precomputed_dispersions,filtering=filtering),
                      series=series,precomputed_dispersions=precomputed_dispersions,filtering=filtering)
    stopCluster(cl)
    rm('cl')
  }
  fits <- bind_rows(fits)
  return(fits)
}


custom_processors <- function(dat,pass){
  # Iowa and Indiana
  # if (unique(dat$administrative_area_level_2)%in%c("Iowa", "Indiana", "Kentucky")){
  
  if (unique(dat$administrative_area_level_1)=="United States"){
    if (!is.na(unique(dat$administrative_area_level_2))){
      if (unique(dat$administrative_area_level_2 != "Washington")){
        pass <- pass & (dat$date > ymd("2020-02-28") )
      }
    }
    # Insert new US custom processors here
  }
  return(pass)
}

outlier_detection <- function(x){
  if (sum(x, na.rm=TRUE) < 100 &  sum(x==0, na.rm=TRUE) > length(x)*0.5){
    return(x)
  } else {
    res <- NULL
    tryCatch({
      res <- tso(ts(x),
                 type="TC", delta=0.1, maxit.iloop = 100, maxit.oloop = 10,  
                 #tsmethod = "auto.arima", args.tsmethod = list(allowdrift = FALSE, ic = "bic", stationary=TRUE),
                 tsmethod="arima", args.tsmethod=list(order=c(1,1,2), method="ML", transform.pars=TRUE, optim.method="BFGS"),
                 cval=4)
    }, error = function(err){
      res <- tryCatch(tso(ts(x),
                          type="TC", delta=0.1, maxit.iloop = 100, maxit.oloop = 10, 
                          #tsmethod = "auto.arima", args.tsmethod = list(allowdrift = FALSE, ic = "bic", stationary=TRUE),
                          tsmethod="arima", args.tsmethod=list(order=c(1,1,2), method="ML", transform.pars=FALSE),
                          cval=4),error=function(e) NULL)
    })
    if (is.null(res)){
      return(x)
    } else {
      res <- res$outliers %>% 
        filter(tstat > 10)
      x[res$ind] <- NA
      return(x)
    }
  }
}

outlier_detection_par <- function(x,id,cl){
  dd <- data.frame('x'=x,'id'=id)
  clusterExport(cl,'dd',envir = environment())
  ids <- unique(id)
  outdet <- function(id,dd.=dd) outlier_detection(dd$x[dd$id==id])
  x <- parSapply(cl,ids,outdet) %>% unlist
  return(x)
}


forecast <- function(nyc,uk){
  nyc_time <- max(nyc$outbreak_time,na.rm=T)
  max_time <- max(uk$outbreak_time,na.rm=T)
  
  ### need to update x from nyc_time+1 to max_time.
  
  ### To do so, we need to project NYC growth rates assuming a parallel trajectory to UK
  ukk <- uk[outbreak_time>nyc_time,c('outbreak_time','growth_rate','p2.5_growth_rate','p97.5_growth_rate')]

  fit <- gam(growth_rate~s(outbreak_time),data=ukk)
  fit2.5 <- gam(p2.5_growth_rate~s(outbreak_time),data=ukk)
  fit97.5 <- gam(p97.5_growth_rate~s(outbreak_time),data=ukk)
  init <- nyc[.N,c('p2.5_growth_rate','growth_rate','p97.5_growth_rate','mean_position')]
  
  last <- data.frame('outbreak_time'=nyc_time)
  dum <- nyc[1:(max_time-nyc_time)]
  dum[,outbreak_time:=(nyc_time+1):max_time]
  dum[,date:=max(nyc$date)+1:.N]
  
  ### predicting NYC growth rates: UK_gam - UK_intercept + NY_latest_value
  
  dum$growth_rate <- fit$fitted.values-c(predict(fit,newdata = last))+init$growth_rate
  dum$p2.5_growth_rate <- fit2.5$fitted.values-c(predict(fit2.5,newdata=last))+init$p2.5_growth_rate
  dum$p97.5_growth_rate <- fit97.5$fitted.values-c(predict(fit97.5,newdata=last))+init$p97.5_growth_rate
  
  dum$mean_position <- init$mean_position+cumsum(dum$growth_rate)
  dum$p2.5_position <- init$mean_position+cumsum(dum$p2.5_growth_rate)
  dum$p97.5_position <- init$mean_position+cumsum(dum$p97.5_growth_rate)
  dum$rm <- NA
  dum$new_confirmed <- NA
  dum$raw_cases <- NA
  dum$forecast <- TRUE
  nyc$forecast <- FALSE
  
  dum <- rbind(nyc[.N],dum)
  dum[1,p2.5_position:=mean_position]
  dum[1,p97.5_position:=mean_position] 
  ## stitching together so the "forecast" includes the last observation and fans out.
  ## This isn't exactly right, since our last observation's estimated mean position
  ## may not be the true position, but we'll assume these errors in mean_position
  ## are small relative to the propagated uncertainty in growth rates
  
  return(rbind(nyc,dum))
}


glue_ma <- function(nyc){
  ix <- which(nyc$forecast)
  x <- nyc$raw_cases
  x[ix] <- exp(nyc[ix,mean_position]) ## Forecast exp(mean_position)
  
  rm <- nyc$rm
  rm[ix] <- round(frollmean(x,7,align='right',na.rm=T))[ix]
  nyc$rm <- rm
  
  rm2.5 <- rm
  x[ix] <- round(exp(nyc[ix,p2.5_position])) ## replace forecasted observations with 
  rm2.5[ix] <- round(frollmean(x,7,align='right',na.rm=T))[ix]
  nyc$p2.5_rm <- rm2.5
  
  
  rm97.5 <- rm
  x[ix] <- round(exp(nyc[ix,p97.5_position]))
  rm97.5[ix] <- round(frollmean(x,7,align='right',na.rm=T))[ix]
  nyc$p97.5_rm <- rm97.5
  
  return(nyc)
}
