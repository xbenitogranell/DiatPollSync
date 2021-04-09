## Functions from Blas Benito (https://github.com/BlasBenito/BaleFire) https://zenodo.org/record/2599859#.YCTzbGhKiUk 


##-----------------------------------------------------------------------------
# Function for testing for autocorrelation of residuals in a model using Moran's I, then fitting an exponential model to the Moran's I estimates. 

# Returns the range (i.e distance (in years) at which autocorrelation has an effect, defined as distance at which exponential curve predicts Moran's I to go below 0.1
# -----------------------------------------------------------------------------

acfTest = function(x, y.model){
  
  a <- NA
  # x = age vector
  # y.model = model
  
  
  par(mfrow=c(2,1))
  plot(x, residuals(y.model), xlab = "Age (cal yr BP)", ylab ="PC1 residuals")
  abline(h=0, lty=2)
  
  y= residuals(y.model)
  
  ncf.cor = correlog(x, y, y, increment=20, quiet=T, resamp=1000, na.rm=T)
  mIcoef = ncf.cor$correlation
  
  # Truncate values greater than 1/ less than -1 (see Duitilleul et al. 2012, pp. 533)	
  mIcoef[mIcoef>1]<-1
  mIcoef[mIcoef<(-1)]<--1
  
  acfResult = mIcoef
  ageDif= ncf.cor$mean.of.class
  
  acfcheckResult= ageDif
  acfP = ncf.cor$p
  
  # Put the acf scores into a matrix- each column is a time window containing a set of Moran's I scores where appropriate.
  
  acfResult1 = acfResult	
  acfResult1[acfP>0.05]<-0
  
  binACF <-c(1, acfResult[1:30])
  binTime<-c(0, ageDif[1:30])
  
  plot(binTime, binACF, pch=20, col="red", ylab="Correlation", xlab= "Distance (years)", xlim=c(0,300))
  points(binTime, binACF, pch=1, col="blue")
  abline(h=0, lty=2)
  
  
  # Fit exponential decay curve where 0 years =1, asymptote = 0. 
  b0=0 
  b1=1
  
  fit.nls<- tryCatch(fit.nls<-nls(binACF~b0+b1*exp(b2*binTime),start=c(b2=-0.02)), error=function(e) a)
  
  b<-class(fit.nls)
  if(b== "nls"){
    
    # summary(fit.nls)$coefficients[1]
    predict.x <- seq(0, 300)
    predict.y<- predict(fit.nls, newdata= list(binTime=predict.x))
    lines(predict.x, predict.y, col="red")
    
    which(predict.y<0.1)[1]
    
  } else {
    
    fit.nls	
  } 
}


# -----------------------------------------------------------------------------
# Function for binning samples in core.
# -----------------------------------------------------------------------------

binFunc = function(xDF, Ages, binWidth, minBin, maxBin){
  
  coreMat = xDF
  
  coreBin = seq(minBin, maxBin, binWidth)
  sampleBin = vector("list", length(coreBin)-1)
  
  for(i in 1:length(sampleBin)){
    toBin = which(Ages>=coreBin[i] & Ages<coreBin[i+1])
    sampleBin[[i]] = coreMat[toBin,]
  }
  
  coreRes = matrix(NA, ncol=ncol(xDF), nrow=length(sampleBin))
  if(ncol(xDF) > 1) {
    for(i in 1:length(sampleBin)) coreRes[i,] = apply(sampleBin[[i]], 2, function(x)mean(x, na.rm=T))
  } else { 
    for(i in 1:length(sampleBin)) coreRes[i,] = mean(sampleBin[[i]], na.rm=T)
  }
  colnames(coreRes) = colnames(xDF)
  rownames(coreRes) = coreBin[-1]
  
  coreRes <- as.data.frame(coreRes)
  coreRes
  
}



#FUNCTION TO SCALE DATA
scaleData=function(data, old.max, old.min, new.max, new.min){
  result=((data - old.min) / (old.max - old.min)) * (new.max - new.min) + new.min
  return(result)
}

#TO CHECK IF NUMBER IS ODD
is.odd = function(x) x %% 2 != 0


#FUNCTION TO MERGE DATASETS AND INTERPOLATE THEM TO A REGULAR GRID
interpolateDatasets<-function(datasets.list, age.column.name, interpolation.time.step){

  #computing age ranges
  age.ranges=sapply(datasets.list, FUN=function(x) range(x[, age.column.name]))
  #min of maximum ages
  min.age=round(max(age.ranges[1,]), 1)
  #max of minimum ages
  max.age=round(min(age.ranges[2,]), 1)

  #subsetting dataframes in list
  datasets.list=lapply(datasets.list, function(x) x[x[, age.column.name] >= min.age & x[, age.column.name] <= max.age, ])

  #reference data
  reference.age = seq(min.age, max.age, by=interpolation.time.step)

  #looping through datasets to interpolate
  for (dataset.to.interpolate in names(datasets.list)){

    #getting the dataset
    temp = datasets.list[[dataset.to.interpolate]]

    #removing age from the colnames list
    colnames.temp = colnames(temp)
    colnames.temp = colnames.temp[which(colnames.temp != age.column.name)]

    #empty dataset to store interpolation
    temp.interpolated = data.frame(age=reference.age)

    #iterating through columns
    for (column.to.interpolate in colnames.temp){

      #do not interpolate non-numeric columns
      if (is.numeric(temp[, column.to.interpolate])==FALSE){
        temp.interpolated[, column.to.interpolate]=temp[, column.to.interpolate]
        next
      }

      #interpolation
      interpolation.formula = as.formula(paste(column.to.interpolate, "~", age.column.name, sep=" "))

      #iteration through span values untill R-squared equals 1
      #span.values=seq(0.01, -0.000001, by=-0.000001) #original code span values should take at least n rows
      span.values=seq(60, -0.650, by=-0.650) #for pollen (&agropastoralism) PrC data
      for(span in span.values){
        interpolation.function = loess(interpolation.formula, data=temp, span=span, control=loess.control(surface="direct"))

        #check fit
        if(cor(interpolation.function$fitted, temp[, column.to.interpolate]) >=  0.99999999){break}

      }

      print(paste("Correlation between observed and interpolated data = ", cor(interpolation.function$fitted, temp[, column.to.interpolate]), sep=""))

      interpolation.function = loess(interpolation.formula, data=temp, span=0.02, control=loess.control(surface="direct"))
      interpolation.result = predict(interpolation.function, newdata=reference.age, se=FALSE)

      #constraining the range of the interpolation result to the range of the reference data
      interpolation.range=range(temp[, column.to.interpolate])
      interpolation.result[interpolation.result < interpolation.range[1]] = interpolation.range[1]
      interpolation.result[interpolation.result > interpolation.range[2]] = interpolation.range[2]

      #putting the interpolated data back in place
      temp.interpolated[, column.to.interpolate]=interpolation.result

    }#end of iteration through columns

    #removing the age column
    temp.interpolated[, age.column.name]=NULL

    #putting the data back in the list
    datasets.list[[dataset.to.interpolate]] = temp.interpolated

  }#end of iterations through datasets

  #same rows?
  nrow.datasets=sapply(datasets.list, FUN=function(x) nrow(x))
  if (length(unique(nrow.datasets))==1){

    #remove age from all dataframes
    datasets.list=lapply(datasets.list, function(x) { x[, age.column.name] = NULL; x })

    #put dataframes together
    output.dataframe = do.call("cbind", datasets.list) #changes names

  } else {
    stop("Resulting datasets don't have the same number of rows, there's something wrong with something.")
  }

  #add reference.age
  output.dataframe = data.frame(age=reference.age, output.dataframe)

  return(output.dataframe)

}


#GENERATES LAGGED CHARCOAL DATA BEFORE EVERY ERICA SAMPLE
backwardLags <- function(lags, reference.data, data.to.lag){
  
  #df to store the lagged data
  lag.data = data.frame(diatom_deriv=double(), pollen_deriv=double(), lag=integer())
  
  #iterates through erica lines
  for (diat.case in 1:nrow(diat)){
    
    #take a line of the diat dataframe and replicate it as many times as lags are
    diat.value = rep(diat[diat.case, "diatom_deriv"], max(lags))

    #get the age of the replicated line
    diat.case.age=diat[diat.case, "age"]
    
    #diat.case.age.plus.lags=round((diat.case.age + (diff(pollen$age)[1] * max(lags))), 2)
    diat.case.age.plus.lags=round((diat.case.age + (diff(agropastolarism$age)[1] * max(lags))), 2)

    #if beyond maximum age
    if (diat.case.age.plus.lags > max(pollen$age)){break}
    #if (diat.case.age.plus.lags > max(agropastolarism$age)){break}
    
    #get from pollen the lines with age > age.diatoms && age <= age.diatoms + lags
    pollen.temp=pollen[which(pollen$age > diat.case.age & pollen$age <= diat.case.age.plus.lags), "pollen_deriv"]
    #pollen.temp=agropastolarism[which(agropastolarism$age > diat.case.age & agropastolarism$age <= diat.case.age.plus.lags), "pollen_deriv"]
    
    #put the data together
    pollen.temp=data.frame(diatom_deriv=diat.value, pollen_deriv=pollen.temp, lag=lags)

    #put them in the final table
    lag.data=rbind(lag.data, pollen.temp)

  }#end of iterations
  
  #remove stuff we don't need
  rm(pollen.temp, diat.case, diat.case.age, diat.case.age.plus.lags, diat.value)

  #order by lag
  lag.data=lag.data[order(lag.data$lag),]
  
  #standardize data
  #get lags column
  lags.column=lag.data$lag
  
  #standardize
  lag.data=scale(lag.data[, c("diatom_deriv", "pollen_deriv")])
  
  #add lag
  lag.data=data.frame(lag.data, lag=lags.column)
  
  #lag as factor
  lag.data$lag=as.factor(lag.data$lag)
  
  return(lag.data)
}



#FUNCTION TO INSERT MINOR TICKS IN GGPLOT
insert_minor <- function(major_labs, n_minor) {labs =
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}


#PLOT LAGS
plotLags <- function(forward.lagged.data, backward.lagged.data, lags, filename){

  pdf(filename, height=6, width=12)
  for (lag in lags){

    #temp data
    temp.future=forward.lagged.data[forward.lagged.data$lag==lag, ]
    temp.past=backward.lagged.data[backward.lagged.data$lag==lag, ]

    plot.past = ggplot(data=temp.past, aes(x=ericaceae.par, y=charcoal.acc.rate)) + geom_point(shape=21, fill="gray50", color="black", size=3, alpha=0.5) +
      xlab("Erica PAR (pollen grains/cm2yr)") +
      ylab("CHAR (particles/cm2 yr)") +
      ggtitle(paste("Erica vs. CHAR, backward lag ", lag, sep="")) +
      theme(text=element_text(size=12), plot.title=element_text(size = 16))

    plot.future = ggplot(data=temp.future, aes(x=ericaceae.par, y=charcoal.acc.rate)) + geom_point(shape=21, fill="gray50", color="black", size=3, alpha=0.5) +
      xlab("Erica PAR (pollen grains/cm2yr)") +
      ylab("CHAR (particles/cm2 yr)") +
      ggtitle(paste("Erica vs. CHAR, forward lag ", lag, sep="")) +
      theme(text=element_text(size=12), plot.title=element_text(size = 16))

    print(cowplot::plot_grid(plot.past, plot.future, align="h", ncol=2))
  }
  dev.off()

}


#COMPUTES GLS MODELS BETWEEN A RESPONSE AND A VARIABLE IN A LAGGED DATASET. PROVIDES A DATAFRAME
modelLagData <- function(model.formula, lagged.data){

  lags=as.numeric(sort(unique(lagged.data$lag)))
  model.formula=as.formula(model.formula)
  response = all.vars(model.formula)[1]

  #list to store results
  results.list = list()

  #fitting a model per lag
  for (lag in lags){
    results.list[[lag]] = gls(model.formula, data=lagged.data[lagged.data$lag==lag,])
  }

  #list to store pseudo R2
  results.list.R2=list()

  #computing pseudo R2 per lag
  for (lag in lags){
    results.list.R2[[lag]] = cor(lagged.data[lagged.data$lag==lag, response], predict(results.list[[lag]]))^2
  }

  #gathering coefficients
  results.list.coef = lapply(results.list, function(x){summary(x)$tTable[2,1]})

  #gathering the standard error of the coefficients
  results.list.coef.se = lapply(results.list, function(x){summary(x)$tTable[2,2]})

  #gathering p-value of coefficients
  results.list.pvalue = lapply(results.list, function(x){summary(x)$tTable[2,4]})

  #to data frame
  output.df = as.data.frame(do.call(rbind, results.list.coef))
  output.df = data.frame(lag=as.numeric(rownames(output.df)), Coefficient=output.df$V1)
  output.df[, "p-value"] = as.data.frame(do.call(rbind, results.list.pvalue))$V1
  output.df$R2 = as.data.frame(do.call(rbind, results.list.R2))$V1
  #output.df$lag=output.df$lag*10
  output.df$lag=output.df$lag*24 #temporal resolution
  #output.df$lag=output.df$lag*60 #temporal resolution
  

  #se of coefficients
  output.df.se = as.data.frame(do.call(rbind, results.list.coef.se))
  output.df.se$lower = output.df$Coefficient - output.df.se$V1
  output.df.se$upper = output.df$Coefficient + output.df.se$V1
  output.df.se$V1 = NULL

  #to long format for plotting
  output.df.long = gather(output.df, variable, value, 2:ncol(output.df))

  #adding the errors
  output.df.long$lower = c(output.df.se$lower, output.df.long[output.df.long$variable=="p-value", "value"], output.df.long[output.df.long$variable=="R2", "value"])
  output.df.long$upper = c(output.df.se$upper, output.df.long[output.df.long$variable=="p-value", "value"], output.df.long[output.df.long$variable=="R2", "value"])

  return(output.df.long)
}



#COMPUTING NULL
#randomizing lag column of both datasets 999 times and computing model with randomized data
modelRandomLagData <- function(lagged.data, model.formula, iterations){

  #getting response column
  response = all.vars(as.formula(model.formula))[1]

  #list to store results
  results = list()

  #computing one model per iteration
  for(i in 1:iterations){

    #randomization of the response column
    lagged.data[, response] = lagged.data[sample(1:nrow(lagged.data)), response]

    #computing model and storing results in list
    results[[i]] = modelLagData(model.formula=model.formula, lagged.data=lagged.data)

  }#end of iterations

  #preparing the data
  #getting lags and variable names
  null.model.df = data.frame(lag=results[[1]]$lag, variable=results[[1]]$variable, stringsAsFactors = FALSE)

  #getting the value column of each dataframe in each list
  null.model.values = sapply(results, `[[`, 3)

  null.model.df$value = apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.5)
  null.model.df$upper = apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.95)
  null.model.df$lower = apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.05)

  return(null.model.df)

}


#FUNCTION TO PLOT MODELING RESULTS
plotModelOutput <- function(backward.results, backward.results.random, width, height, title.size, text.size, filename){

  #axes limits
  max.lag = max(c(backward.results$lag))
  max.coefficient = round(max(c(backward.results[backward.results$variable=="Coefficient", "value"], backward.results.random[backward.results.random$variable=="Coefficient", "upper"])) + 0.1, 1)
  min.coefficient = round(min(c(backward.results[backward.results$variable=="Coefficient", "value"], backward.results.random[backward.results.random$variable=="Coefficient", "lower"])) - 0.1, 1)
  max.R2 = round(max(c(backward.results[backward.results$variable=="R2", "value"])), 1)

  #separating pvalues
  # forward.results.pvalue = forward.results[forward.results$variable=="p-value", c("lag", "value")]
  # backward.results.pvalue = backward.results[backward.results$variable=="p-value", c("lag", "value")]

  #reference color scale
  library(viridis)
  viridis.colors = viridis(10, option="D")

  #BACKWARD PLOT
  backward.plot.coefficient = ggplot(data=subset(backward.results, variable=="Coefficient"), aes(x=lag, y=value)) +
    geom_ribbon(data=subset(backward.results.random, variable=="Coefficient"), aes(ymin=lower, ymax=upper), alpha=0.3, fill=viridis.colors[10]) +
    geom_line(data=subset(backward.results.random, variable=="Coefficient"), aes(x=lag, y=value), alpha=0.6, color=viridis.colors[10], size=1) +
    geom_hline(yintercept=0, color="black", linetype=3) +
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.5, fill=viridis.colors[5]) +
    geom_line(size=1.5, color=viridis.colors[1]) +
    ggtitle(expression("Pollen" %->% "Diatoms")) +
    #ggtitle(expression("Agropastoralism indicators" %->% "Diatoms")) +
    theme(legend.position="none") +
    xlab("") +
    ylab("Standardized coefficient") +
    scale_y_continuous(breaks=seq(min.coefficient, max.coefficient, by=0.2)) +
    #scale_x_reverse(limits=c(max.lag, 10), breaks=c(10, seq(100, max.lag, by=100))) +
    scale_x_reverse(limits=c(max.lag, diff(pollen$age)[1]), 
                    breaks=c(diff(pollen$age)[1], seq(100,max.lag, by=960))) +
    #scale_x_reverse(limits=c(max.lag, diff(agropastolarism$age)[1]), breaks=c(diff(agropastolarism$age)[1], seq(40, max.lag, by=diff(agropastolarism$age)[1]))) +
    theme(axis.text = element_text(size=text.size),
          plot.title = element_text(size = title.size),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x = element_blank()) +
    coord_cartesian(ylim = c(min.coefficient, max.coefficient + 0.1))

  
  backward.plot.R2 = ggplot(data=subset(backward.results, variable=="R2"), aes(x=lag, y=value, group=variable)) +
    geom_ribbon(data=subset(backward.results.random, variable=="R2"), aes(ymin=lower, ymax=upper), alpha=0.3, fill=viridis.colors[10]) +
    geom_line(data=subset(backward.results.random, variable=="R2"), aes(x=lag, y=value), alpha=0.6, color=viridis.colors[10], size=1.5) +
    geom_line(size=1.5, color=viridis.colors[1]) +
    theme(legend.position="none") +
    xlab("Years (before Diatom samples)") +
    ylab("Pseudo R squared") +
    scale_y_continuous(breaks=seq(0, max.R2, by=0.1)) +
    #scale_x_reverse(limits=c(max.lag, 10), breaks=c(10, seq(100, max.lag, by=100))) +
    #scale_x_reverse(limits=c(max.lag, diff(agropastolarism$age)[1]), breaks=c(diff(agropastolarism$age)[1], seq(40, max.lag, by=diff(agropastolarism$age)[1]))) +
    scale_x_reverse(limits=c(max.lag, diff(pollen$age)[1]), 
                    breaks=c(diff(pollen$age)[1], seq(100,max.lag, by=960))) +
    theme(axis.text = element_text(size=text.size),
          axis.text.x = element_text(size=text.size),
          plot.title = element_text(size = title.size),
          plot.margin = unit(c(0.2, 0.5, 0, 0), "cm")) +
    coord_cartesian(ylim = c(0, max.R2 + 0.05))

  
  first_col = plot_grid(backward.plot.coefficient, backward.plot.R2, ncol = 1, rel_heights = c(1, 1), align="v") + theme(plot.margin = unit(c(0.5, -1, 0.5, 0.5), "cm"))

  #complete.plot = plot_grid(first_col, NULL, second_col, ncol = 3, rel_widths = c(1, 0, 1), align="h")

  #print(complete.plot)
  print(first_col)
  
  ggsave(filename = filename, width=width, height=height)

}




