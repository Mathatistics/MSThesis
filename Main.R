
## ----frontMatter, child="frontMatter.Rnw"--------------------------------


## ----LoadingPkgs, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
req.package<-c("MASS", "car", "pls", "xtable", "grid", "gridExtra", "knitr", "leaps", "zoo", "gdata","ridge", "plyr", "dplyr", "ggplot2", "reshape2", "scales","mixlm")
lapply(req.package, require, character.only=TRUE, quietly = T, warn.conflicts = F)




## ----setup, include=FALSE, cache=FALSE, echo=TRUE------------------------
opts_chunk$set(fig.path='Include', fig.align='center')
render_listings()
setwd('~/Dropbox/UMB/Thesis/MSThesis/')
Sys.setenv(TEXINPUTS=getwd(),
           BIBINPUTS=getwd(),
           BSTINPUTS=getwd())
#data.path<-path.expand(file.path(dirname(dirname(getwd())), "Datasets", "CompleteDataSet.xlsx"))
data.path<-path.expand(file.path(dirname(getwd()), "Datasets", "CompleteDataSet.xlsx"))


## ----readFun, child="Include/functions.Rnw"------------------------------


## ----functions, echo=FALSE, cache=FALSE, warning=FALSE-------------------

## Setting up Crisis Period
cp.cat<-function(dateVec){
cp.col<-ifelse(dateVec<cperiod[1] | dateVec>cperiod[2],
               "Normal Period", 
               "Crisis Period")
return(cp.col)
}

## Timeseries plot
plotTS<-function(dataSet, dateVarColIdx, nc){
  plt<-ggplot(melt(dataSet, dateVarColIdx), aes(Date, (value/100)))
  plt<-plt+geom_line()
  plt<-plt+facet_wrap(~variable, 
                      ncol=nc, 
                      scale="free_y")
  plt<-plt+theme_bw()
  plt<-plt+theme(text=element_text(size=12))
  plt<-plt+labs(x="Date (Monthly)", y="Value (NOK hundreds)")
  return(plt)
}

## Plotting Model Coefficients with their state of significance
test.plot<-function(model, alpha=0.05){
  .e<-environment()
  coef.matrix<-data.frame(summary(model)$coef)
  names(coef.matrix)<-c("Estimate", "StdError", "t.value", "p.value")
  idx<-order(row.names(coef.matrix))
  cp<-ggplot(coef.matrix[idx,], aes(x=row.names(coef.matrix[idx,]), y=t.value), environment = .e)
  cp<-cp+geom_bar(stat="identity", position = "identity",
                  fill=ifelse(coef.matrix[idx,"p.value"]<alpha, "coral3", "cornflowerblue"))
  cp<-cp+geom_text(aes(y=ifelse(coef.matrix[idx, "t.value"]>0,t.value+0.7, t.value-0.7), 
                       label=round(coef.matrix[idx,"Estimate"], 2)), angle=45, size=5)
  cp<-cp+theme_bw()+labs(x="", y="T-Value")
  cp<-cp+theme(axis.text.x=element_text(angle=90, hjust=1))
  cp<-cp+theme(text=element_text(size=20))
  cp<-cp+scale_fill_manual("Status", values=c("firebrick2", "dodgerblue3"), 
                           labels=c("Significant", "Non-Significant"))
  cp<-cp+geom_hline(yintercept=c(-1,1)*qt(alpha/2, df = abs(diff(dim(model$model[,-1]))), lower.tail = F), 
                    color="red", linetype="dashed")
  cp<-cp+theme(legend.title=element_blank(), 
               legend.position=c(0.8, 0.2))
  cp<-cp+geom_hline(yintercept=0, color="black", size=.2)
  return(cp)
}

## Fitting Linear Model
fit.model<-function(Model, yVar, xVars, dataSet, scaling=TRUE){
  model<-match.fun(Model)
  formula<-as.formula(paste(yVar, paste(xVars, collapse="+"), sep="~"))
  if(scaling){
      model<-model(formula, data=dataSet, scale=TRUE)
  }else{
      model<-model(formula, data=dataSet)
  }
  return(list(formula=formula, model=model, dataset=dataSet))
}


## Diagnostic Plot using GGPlot
diagPlot<-function(model, cp.color){
  p1<-ggplot(model, aes(.fitted, .resid))+geom_point(aes_string(color=cp.color))
  p1<-p1+stat_smooth(method="loess")
  p1<-p1+geom_hline(yintercept=0, col="red", linetype="dashed")
  p1<-p1+xlab("Fitted values")+ylab("Residuals")
  p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_bw()
  
  ## qline slope and intercept
  qline<-ldply(data.frame(res=stdres(mdl.ft$linear$model)), function(x){
      slope = (quantile(x,p=.75)-quantile(x,.25))/(qnorm(.75)-qnorm(.25))
      intercept = quantile(x,.25) - slope*qnorm(.25)
      data.frame(slope, intercept)})
  
  p2<-ggplot(model, aes(sample=.stdresid))+stat_qq(aes_string(color=cp.color))
  p2<-p2+geom_abline(data = qline, aes(slope, intercept))+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
  p2<-p2+ggtitle("Normal Q-Q")+theme_bw()
  
  p3<-ggplot(model, aes(.fitted, sqrt(abs(.stdresid))))+geom_point(na.rm=TRUE, aes_string(color=cp.color))
  p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
  p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
  p3<-p3+ggtitle("Scale-Location")+theme_bw()
  
  p4<-ggplot(model, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity", aes_string(fill=cp.color))
  p4<-p4+xlab("Obs. Number")+ylab("Cook's distance")
  p4<-p4+geom_text(aes(x=which.max(.cooksd), 
                   y = max(.cooksd), 
                   label=format(baseTable[which.max(.cooksd), "Date"], "%b %Y")),
                   size=4)
  p4<-p4+ggtitle("Cook's distance")+theme_bw()
  
  p5<-ggplot(model, aes(.hat, .stdresid))
  p5<-p5+geom_point(aes_string(color=cp.color, size=".cooksd"), na.rm=TRUE)
  p5<-p5+stat_smooth(method="loess", na.rm=TRUE)
  p5<-p5+xlab("Leverage")+ylab("Standardized Residuals")
  p5<-p5+ggtitle("Residual vs Leverage Plot")
  p5<-p5+scale_size_continuous("Cook's Distance", range=c(1,5))
  p5<-p5+theme_bw()+theme(legend.position="bottom")
  
  p6<-ggplot(model, aes(.hat, .cooksd))+geom_point(na.rm=TRUE, aes_string(color=cp.color))+stat_smooth(method="loess", na.rm=TRUE)
  p6<-p6+xlab("Leverage hii")+ylab("Cook's Distance")
  p6<-p6+ggtitle("Cook's dist vs Leverage hii/(1-hii)")
  p6<-p6+geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")
  p6<-p6+theme_bw()
  
  return(list(rvfPlot=p1, qqPlot=p2, sclLocPlot=p3, cdPlot=p4, rvlevPlot=p5, cvlPlot=p6))
}

## Generate summary plot from a fitted model to annotate other plot
sumryBlock<-function(model){
  return(paste("R-Sq = ",signif(summary(model)$r.squared, 3),
               "\nAdj R-Sq =",signif(summary(model)$adj.r.squared, 3),
               "\nSigma =",signif(summary(model)$sigma, 3),
               "\nF =",signif(as.vector(summary(model)$fstatistic[1]), 4),
               paste("(",paste(as.vector(summary(mdl.ft$cp.model$model)$f[2:3]), collapse=','),")", sep="")
         ))
}

model.sumry<-function(model, call=TRUE, coefMat=TRUE, sumry=TRUE){
    if(!"lm"%in%class(model)){
        stop("Model should be of class 'lm'.\n")
    }
    else{
        s<-summary(model)$sigma
        df<-summary(model)$df
        r.sq<-summary(model)$r.squared
        adj.r.sq<-summary(model)$adj.r.squared
        f<-summary(model)$fstatistic[1]
        f.df.num<-summary(model)$fstatistic[2]
        f.df.den<-summary(model)$fstatistic[3]
        if(call){
            print(summary(model)$call)
            cat("\n")
        }
        if(coefMat){
            printCoefmat(summary(model)$coef, digits = 3)
        }
        if(sumry){
            data.frame(Sigma=summary(model)$sigma, 
                       R.Sq=summary(model)$r.squared, 
                       R.Sq.adj=summary(model)$adj.r.squared, 
                       F.value=summary(model)$fstatistic[1], 
                       df=paste(summary(model)$fstatistic[2:3], collapse=","), 
                       p.value=pf(summary(model)$fstatistic[1], 
                                  summary(model)$fstatistic[2],
                                  summary(model)$fstatistic[3], 
                                  lower.tail = FALSE))
        }
    }
}

vifPlot<-function(model){
    if("lm"%nin%class(model)){
        stop("Model should be of class 'lm'.")
    }else{
        coef<-names(vif(model))
        vif<-as.vector(vif(model))
        mdl.label<-ifelse(label(model)=="", deparse(substitute(model)), label(model))
        vifMat<-data.frame(coef, vif)
        p<-ggplot(vifMat, aes(coef, vif))
        p<-p+geom_bar(stat="identity", color="black", fill=NA)+theme_bw()
        p<-p+ggtitle(label = paste("Variance Inflation Function plot\nModel:", mdl.label))
        if(length(coef)>5){
            p<-p+theme(axis.text.x=element_text(hjust=1, angle=90))
        }
        return(p)
    }
}

addline_format <- function(x,...){
    gsub('\\s','\n',x)
}


## Function to perform cross-validation splitting into 12 consecutive segments on Linear model and its subsets
makeFormula<-function(x.var, y.var){
    formula<-paste(y.var, paste(x.var, collapse="+"), sep="~")
    return(formula)
}
mdl.cv<-function(dataSet, x.var, y.var, model="lm", step=FALSE, criteria=NULL, split=12, lmd=NULL){
    segment<-split(1:nrow(dataSet), ceiling(1:nrow(dataSet)/split))
    formula=makeFormula(x.var, y.var)
    mdl<-list()
    predVec<-rep(NA, nrow(dataSet))
    errVec<-rep(NA, nrow(dataSet))
    
    for(i in seq_along(segment)){
        dataset<-dataSet[-segment[[i]],]
        testset<-dataSet[segment[[i]],]
        if(step & model=="lm"){
            if(!criteria %in% c("AIC", "BIC", "Cp", "R2adj", "forward", "backward")){
                stop("Please! enter the correct criteria")
            }else{
                require(leaps)
                if(criteria=="Cp"){
                    ## Model selected by Mallows Cp Criteria
                    cp.leaps<-leaps(x=dataset[,x.var],
                                    y=dataset[,y.var],
                                    method="Cp", nbest = 1, names = x.var)
                    # Model fitting
                    cp.which<-names(which(cp.leaps$which[which.min(cp.leaps$Cp),]))
                    formula<-makeFormula(cp.which, y.var)
                    mdl[[i]]<-lm(formula, data=dataset)
                }else if(criteria=="R2adj"){
                    ## Model selected by R2adj Criteria
                    r2adj.leaps<-leaps(x=dataset[,x.var],
                                       y=dataset[,y.var],
                                       method="adjr2", nbest = 1, names=x.var)
                   # Model fitting
                    r2.which<-names(which(r2adj.leaps$which[which.max(r2adj.leaps$adjr2),]))
                    formula<-makeFormula(r2.which, y.var)
                    mdl[[i]]<-lm(formula, data=dataset)
                }else if(criteria=="AIC" | criteria=="BIC"){
                    lmBstSetSmry <- summary(regsubsets(dataset[,x.var],
                                                       dataset[,y.var], 
                                                       nbest = 1, nvmax = length(x.var)))
                    nvars<-apply(lmBstSetSmry$which, 1, sum)
                    bic.vec<-lmBstSetSmry$bic
                    aic.vec<-bic.vec-nvars*log(sum(train))+nvars
                    
                    ## Fitting selected linear model
                    aic.which<-names(which(lmBstSetSmry$which[which.min(aic.vec),]))[-1]
                    bic.which<-names(which(lmBstSetSmry$which[which.min(bic.vec),]))[-1]
                    if(criteria=="AIC"){
                        formula<-makeFormula(aic.which, y.var)
                        mdl[[i]]<-lm(formula, data=dataset)
                    }else if(criteria=="BIC"){
                        formula<-makeFormula(bic.which, y.var)
                        mdl[[i]]<-lm(formula, data=dataset)
                    }
                }else if(criteria=="forward"){
                        require(mixlm)
                        fm.log<-capture.output({
                        mdl[[i]]<- forward(do.call(lm, list(formula, dataset)), alpha = 0.05, full = FALSE)
                    })
                }else if(criteria=="backward"){
                        require(mixlm)
                        fm.log<-capture.output({
                        mdl[[i]]<- backward(do.call(lm, list(formula, dataset)), alpha = 0.05, full = FALSE)
                    })
                }
            }
        }else if(step & model!='lm'){
            stop("Stepwise can only be performed using Linear Model, Please input 'lm' in the model.")
        }else if(model=='lm'){
            mdl[[i]]<-lm(formula, dataset)
        }else if(model=='ridge'){
            require(ridge)
            mdl[[i]]<- linearRidge(formula, dataset, lambda = lmd)
        }else{
            stop("Model can take 'lm' or 'ridge' value.")
        }
        predVec[segment[[i]]]<-predict(mdl[[i]], newdata=testset[,x.var])
        errVec[segment[[i]]]<-testset[,y.var]-predVec[segment[[i]]]
    }
    rmse.cv<-sqrt(1/nrow(dataSet)*sum(errVec^2))
    r2pred<-1-sum(errVec^2)/sum((predVec-mean(dataSet[,y.var]))^2)
    invisible(list(Model=mdl, Predicted=predVec, Error=errVec, rmsep=rmse.cv, r2pred=r2pred))
}

## Grid Arrange with common Legend
grid_arrange_shared_legend <- function(plotList, ncol=2, main=NULL, ...) {
    plots <- plotList
    g <- ggplotGrob(plots[[1]] + 
                        theme(legend.position="bottom", 
                              legend.title=element_blank()))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    plt.lst<-lapply(plots, function(x){
        x + theme(legend.position="none")
    })
    plt.lst$ncol<-ncol
    plt.lst$main<-main
    grid.arrange(
        do.call(arrangeGrob, plt.lst),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}





## ----data-prep, child="Include/DataPreperation.Rnw"----------------------


## ----dataSetup, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
baseTable<-read.xls(data.path, sheet = "FinalData")
baseTable[,1]<-as.Date(baseTable[,1], format="%d/%m/%Y")
baseTable[,"Testrain"]<-as.logical(baseTable[,"Testrain"])
# baseTable1<-baseTable

## Log Transform some variable using log1p() Function
## baseTable[, "ImpOldShip"]<-log1p(baseTable[, "ImpOldShip"])
# baseTable[, "ExpOilPlat"]<-log1p(baseTable[, "ExpOilPlat"])
# baseTable[, "ExpExShipOilPlat"]<-log1p(baseTable[, "ExpExShipOilPlat"])


## Label Variables in baseTable
labelTable<-read.xls(data.path, sheet = "FinalCodeBook", stringsAsFactors=FALSE)
for(i in 1:ncol(baseTable)){
    Hmisc::label(baseTable[,i])<-labelTable[i,2]
    class(baseTable[,i])<-rev(class(baseTable[,i]))
}

# Variable Declaration
y.var<-grep("PerEURO", names(baseTable), value=TRUE)
fin.var<-grep("^CPI|Int", names(baseTable), value=TRUE)
price.var<-grep("^Oil", names(baseTable), value=TRUE)
import.var<-grep("^Imp", names(baseTable), value=TRUE)
export.var<-grep("^Exp", names(baseTable), value=TRUE)
tradeBal.var<-grep("^Tr", names(baseTable), value=TRUE)
expct.var<-grep("^l", names(baseTable), value=TRUE)
y2.var<-grep("ExcCh", names(baseTable), value=TRUE)
season<-grep("season", names(baseTable), value=TRUE)
train<-grep("Testrain", names(baseTable), value=TRUE)

x.var<-c(fin.var, price.var, import.var, export.var, tradeBal.var, expct.var)
# baseTable$Testrain<-baseTable$Date<"2013-01-01"
train<-baseTable[,"Testrain"]

balTot<-balTot<-read.xls(file.path(dirname(data.path), "Balance of Payment Quarterly Data.xlsx"), sheet = "BalTot")
balTot<-balTot[-nrow(balTot),]
balTot$Date<-as.yearqtr(gsub("K", "Q", balTot$Date))

## Crisis Period
cperiod<-c("2007-06-01", "2009-06-01") ## Three Years of crisis Period




## ----chapter1-include, child="Include/Chapter-1.Rnw", eval=TRUE----------




## ----chapter2-include, child="Include/Chapter-2.Rnw", eval=TRUE----------












## ----tsPlotExp, echo=FALSE, fig.height=5, fig.cap="Time Series plot of major exports of Norway", warning=FALSE, error=FALSE----
plotTS(baseTable[,c("Date", ls(baseTable, pattern = "Exp"))], 1, nc=2)










## ----chapter3-include, child="Include/Chapter-3.Rnw", eval=TRUE----------


## ----mdlFitCriteriaPlot, child="mdlFitCriteriaPlot.Rnw"------------------






## ----chapter4-include, child="Include/Chapter-4.Rnw", eval=TRUE----------


## ----sumryTablSetup, echo=FALSE, results='hide'--------------------------
sumryTabl<-t(sapply(baseTable[,c(y.var, x.var)], 
                           function(x){c(min=min(x), 
                                         median=median(x), 
                                         max=max(x), 
                                         mean=mean(x), 
                                         stdev=sd(x))}))
sumryXtable<-xtable(sumryTabl)

## Repeat Table Header Row for longtable ########
addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(0)
addtorow$command  <- c(paste("\\hline \n",
                             "\\endhead \n",
                             "\\hline \n",
                             "{\\footnotesize Continued on next page} \n",
                             "\\endfoot \n",
                             "\\endlastfoot \n",sep=""))
## ------------------------ #########

caption(sumryXtable)<-"Summary Report of all the variables used in this report"
label(sumryXtable)<-"tbl:sumryTabl"







## ----commons, child="Include/Commons.Rnw"--------------------------------


## ----modelFitting, echo=FALSE, results='hide'----------------------------
pls.options(plsralg="oscorespls")
mdl<-c("lm", "pcr", "plsr", "linearRidge")
mdl.ft<-lapply(seq_along(mdl), 
               function(x){
                   do.call(fit.model, list(
                        mdl[x], 
                        y.var, 
                        x.var, 
                        baseTable[train,], 
                        scaling=c(mdl %in% c("plsr", "pcr"))[x]
                   ))
                   })
names(mdl.ft)<-c("linear", "PCR", "PLS", "ridge")

## --------------------------------------------------------------------|
## Model selected by Mallows Cp Criteria
cp.leaps<-leaps(x=mdl.ft$linear$dataset[,x.var],
                y=mdl.ft$linear$dataset[,y.var],
                method="Cp", nbest = 1, names = x.var)

# Prepare for plot
cpdf<-data.frame(p=cp.leaps$size, cp=cp.leaps$Cp)

# Model fitting
cp.which<-names(which(cp.leaps$which[which.min(cp.leaps$Cp),]))
mdl.ft$cp.model<-do.call(fit.model, list("lm", y.var, cp.which, baseTable[train,], scaling = FALSE))

## --------------------------------------------------------------------|
## Model selected by R-sq Adjusted Criteria
r2adj.leaps<-leaps(x=mdl.ft$linear$dataset[,x.var],
                y=mdl.ft$linear$dataset[,y.var],
                method="adjr2", nbest = 1, names=x.var)
# Prepare for plot
r2df<-data.frame(p=r2adj.leaps$size, r2adj=r2adj.leaps$adjr2)

# Model fitting
r2.which<-names(which(r2adj.leaps$which[which.max(r2adj.leaps$adjr2),]))
mdl.ft$r2.model<-do.call(fit.model, list("lm", y.var, r2.which, baseTable[train,], scaling=FALSE))

## --------------------------------------------------------------------|
## Model selected by AIC and BIC criteria
lmBstSetSmry <- summary(regsubsets(mdl.ft$linear$dataset[,x.var],
                                   mdl.ft$linear$dataset[,y.var], 
                                   nbest = 1, nvmax = length(x.var)))
nvars<-apply(lmBstSetSmry$which, 1, sum)
bic.vec<-lmBstSetSmry$bic
aic.vec<-bic.vec-nvars*log(sum(train))+nvars
infoMat<-data.frame(p=nvars, aic=aic.vec, bic=bic.vec)

## Fitting selected linear model
aic.which<-names(which(lmBstSetSmry$which[which.min(aic.vec),]))[-1]
bic.which<-names(which(lmBstSetSmry$which[which.min(bic.vec),]))[-1]

mdl.ft$aicMdl<- do.call(fit.model, list("lm", y.var, aic.which, dataSet = baseTable[train,], scaling = F))
mdl.ft$bicMdl<- do.call(fit.model, list("lm", y.var, bic.which, dataSet = baseTable[train,], scaling = F))

## --------------------------------------------------------------------|
## Forward Selection Model (criteria: level of significance)
fw.model.log <- capture.output(fw.model<-forward(lm(formula = mdl.ft$linear$formula, data=mdl.ft$linear$data), alpha = 0.1, full = FALSE))
mdl.ft$forward<-list(formula=mdl.ft$linear$formula, model=fw.model, data=mdl.ft$linear$data)

## Backward Elimination Model (criteria: level of significance)
bw.model.log<-capture.output(bw.model<-backward(lm(formula = mdl.ft$linear$formula, data=mdl.ft$linear$data), alpha = 0.1, full = FALSE, hierarchy = TRUE))
mdl.ft$backward<-list(formula=mdl.ft$linear$formula, model=bw.model, data=mdl.ft$linear$data)

## --------------------------------------------------------------------|
## Labeling the models
mdl.labels<-c("Linear Model", "Principal Component Regression", "Partial Least Square Regression", "Ridge Regression", "Subset Model (criteria:Mallows Cp)", "Subset Model (criteria:R-sq adjusted)", "Model selected (criteria:AIC)","Model selected (criteria:BIC)", "Forward Selection Model(criteria:F-test)", "Backward Elimination Model (criteria: F-test)")
mdl.prnt.lab<-c("Linear Model", "Principal Component \\\\ Regression", "Partial Least Square \\\\ Regression", "Ridge Regression", "Subset Model \\\\ (criteria:Mallows Cp)", "Subset Model \\\\ (criteria:R-sq adjusted)", "Model selected \\\\ (criteria:AIC)","Model selected (criteria:BIC)", "Forward Selection Model \\\\ (criteria:F-test)", "Backward Elimination Model \\\\ (criteria: F-test)")

for(i in 1:length(mdl.ft)){
    # Label the model
    Hmisc::label(mdl.ft[[i]][[2]])<-mdl.labels[i]
    # Reverse the class
    class(mdl.ft[[i]][[2]])<-rev(class(mdl.ft[[i]][[2]]))
}

## --------------------------------------------------------------------|
## Principal Component Analysis
pc.a<-princomp(baseTable[, x.var], cor = TRUE, )

## --------------------------------------------------------------------|
## Setting up Ridge Parameter lambda
lmd.seq<-seq(0,0.025,0.001)
tuningRidge<-ldply(lmd.seq, function(x){
    rdg.rmsep<-mdl.cv(baseTable[train,], x.var, y.var, 
                      model="ridge", split=12, lmd = x)$rmsep
    rdg.r2pred<-mdl.cv(baseTable[train,], x.var, y.var, 
                       model="ridge", split=12, lmd = x)$r2pred
    data.frame(lmd=x, rmsep=rdg.rmsep, r2pred=rdg.r2pred)
})
tuningRidge<-data.frame(tuningRidge)
lmd<-lmd.seq[which.min(tuningRidge$rmsep)]

## --------------------------------------------------------------------|
## Updating Linear Ridge model with new paramter lmd
mdl.ft$ridge$model<-linearRidge(mdl.ft$ridge$formula, 
                                data=mdl.ft$ridge$dataset, 
                                lambda = lmd)

## Color for crisis period
cperiod.col<-cp.cat(cperiod)



## ----chapter4a-include, child="Include/Chapter-4a.Rnw", eval=TRUE--------


## ----sigCoef, echo=FALSE, warning=FALSE, results='hide'------------------
coefMat<-as.data.frame(summary(mdl.ft$linear$model)$coefficients)
sigVarIdx<-which(coefMat$`Pr(>|t|)`<=0.05)

















## ----chapter4b-include, child="Include/Chapter-4b.Rnw", eval=TRUE--------


## ----pcaSumrySetup, echo=FALSE, results='hide'---------------------------
stdev<-pc.a$sdev
varprop<-pc.a$sdev^2/sum(pc.a$sdev^2)
pcaSumry<-data.frame(cbind( `Comp`=1:length(varprop),
                            `Std.Dev`=stdev, 
                            `Var.Prop`=varprop, 
                            `Cum.Var.Prop`=cumsum(varprop)))
pcaSumry$Comp<-1:nrow(pcaSumry)
pcaSumry1<-xtable(cbind(pcaSumry[1:7,], pcaSumry[8:14,]), digits = 3)
caption(pcaSumry1)<- "Dispersion of data explained by principal components"
label(pcaSumry1)<- "tbl:pcaSumry"
align(pcaSumry1)<- "lrrrr|rrrr"



## ----pcrSumrySetup, echo=FALSE, results='hide'---------------------------
pcr.expVar.x<-cumsum(explvar(mdl.ft$PCR$model))
pcr.expVar.y<-apply(fitted(mdl.ft$PCR$model), 3, var)/var(mdl.ft$PCR$dataset[,y.var])*100
pcrSumry<-data.frame(Comp=1:length(pcr.expVar.x), 
                     X=pcr.expVar.x, 
                     PerEURO=pcr.expVar.y, 
                     row.names = NULL)





## ----chapter4c-include, child="Include/Chapter-4c.Rnw", eval=TRUE--------


## ----plsSumry, echo=FALSE, results='hide'--------------------------------
pls.expVar.x<- cumsum(explvar(mdl.ft$PLS$model))
pls.expVar.y<-apply(fitted(mdl.ft$PLS$model), 3, var)/var(mdl.ft$PCR$dataset[,y.var])*100
plsSumry<-data.frame(Comp=1:length(pls.expVar.x), X=pls.expVar.x, PerEURO=pls.expVar.y, row.names = NULL)



## ----PLSnPCRcomp, echo=FALSE, results='hide'-----------------------------
PLSnPCRcomp<-melt(list(`PCR Model`=list(`Predictor Variable`=pcr.expVar.x, 
                                `Response Variable`=pcr.expVar.y),
                       `PLS Model`=list(`Predictor Variable`=pls.expVar.x, 
                                `Response Variable`=pls.expVar.y)))
names(PLSnPCRcomp)<-c("Variance Explained", "type", "model")
PLSnPCRcomp$Components<-factor(1:length(pcr.expVar.x), levels = 1:length(pcr.expVar.x))






## ----chapter4d-include, child="Include/Chapter-4d.Rnw", eval=TRUE--------


## ----rmsepPLSnPCR, echo=FALSE--------------------------------------------
## Fitting PCR and PLS using Cross-validation
pcr.cv<-pcr(mdl.ft$PCR$formula, data=mdl.ft$PCR$dataset, 
            scale=TRUE, validation="CV", segments=12, 
            segments.type="consecutive")
pls.cv<-plsr(mdl.ft$PCR$formula, data=mdl.ft$PCR$dataset, 
            scale=TRUE, validation="CV", segments=12, 
            segments.type="consecutive")
## RMSEP using Cross-validation
rmsep.pcr<-data.frame(comp=RMSEP(pcr.cv)$comps, 
                      r2pred=as.vector(R2(pcr.cv)$val), 
                      t(sapply(RMSEP(pcr.cv)$comps, 
                               function(x){RMSEP(pcr.cv)$val[,,x+1]})))
rmsep.pls<-data.frame(comp=RMSEP(pls.cv)$comps, 
                      r2pred=as.vector(R2(pls.cv)$val), 
                      t(sapply(RMSEP(pls.cv)$comps, 
                               function(x){RMSEP(pls.cv)$val[,,x+1]})))
rmsep.mat<-melt(list(PCR=rmsep.pcr, PLS=rmsep.pls), 1)

## ----cvStat, echo=FALSE--------------------------------------------------
pcr.sc<-15:17
pls.sc<-6:9

lm.cv<-mdl.cv(baseTable[train,], x.var, y.var)
aic.cv<-mdl.cv(baseTable[train,], x.var, y.var, step = TRUE, criteria = "AIC", split = 12)
bic.cv<-mdl.cv(baseTable[train,], x.var, y.var, step = TRUE, criteria = "BIC", split = 12)
backward.cv<-mdl.cv(baseTable[train,], x.var, y.var, step = TRUE, criteria = "backward", split = 12)
ridge.cv<-mdl.cv(baseTable[train,], x.var, y.var, step=FALSE, split=12, model = "ridge", lmd = lmd)

rmse.cv<-data.frame(RMSEP=c(Linear=lm.cv$rmsep, 
            AICModel=aic.cv$rmsep, 
            BICModel=bic.cv$rmsep,
            BackModel=backward.cv$rmsep,
            Ridge=ridge.cv$rmsep,
            PCR=rmsep.pcr[rmsep.pcr$comp%in%pcr.sc, "adjCV"],
            PLS=rmsep.pls[rmsep.pls$comp%in%pls.sc, "adjCV"]))
r2pred.cv<-data.frame(R2pred=c(Linear=lm.cv$r2pred, 
            AICModel=aic.cv$r2pred, 
            BICModel=bic.cv$r2pred,
            BackModel=backward.cv$r2pred,
            Ridge=ridge.cv$r2pred,
            PCR=rmsep.pcr[rmsep.pcr$comp%in%pcr.sc, "r2pred"],
            PLS=rmsep.pls[rmsep.pls$comp%in%pls.sc, "r2pred"]))
cvStat<-data.frame(rmse.cv, r2pred.cv)
rownames(cvStat)[grep("PCR", rownames(cvStat))]<-paste("PCR.Comp", pcr.sc, sep="")
rownames(cvStat)[grep("PLS", rownames(cvStat))]<-paste("PLS.Comp", pls.sc, sep="")

pls.min.comp<-as.numeric(summarize(cvStat[grep("PLS", rownames(cvStat)), ], pls.sc[which.min(RMSEP)]))
pcr.min.comp<-as.numeric(summarize(cvStat[grep("PCR", rownames(cvStat)), ], pcr.sc[which.min(RMSEP)]))




## ----predMat, echo=FALSE-------------------------------------------------
lm.pred<-predict(mdl.ft$linear$model, 
                      newdata = baseTable[!train, x.var])
pcr.pred<-list()
pls.pred<-list()
pcr.pred<-lapply(pcr.sc, function(x){as.vector(predict(mdl.ft$PCR$model, 
                                   newdata = baseTable[!train, x.var], 
                                   ncomp = x))})
pls.pred<-lapply(pls.sc, function(x){as.vector(predict(mdl.ft$PLS$model,
                                   newdata=baseTable[!train, x.var],
                                   ncomp=x))})
names(pcr.pred)<-paste("Comp",pcr.sc, sep="")
names(pls.pred)<-paste("Comp",pls.sc, sep="")

ridge.pred<-predict(mdl.ft$ridge$model,
                         newdata = baseTable[!train, x.var])
cp.model.pred<-predict(mdl.ft$cp.model$model,
                       newdata=baseTable[!train, x.var])
aicMdl.pred<-predict(mdl.ft$aicMdl$model,
                       newdata=baseTable[!train, x.var])
bicMdl.pred<-predict(mdl.ft$bicMdl$model,
                       newdata=baseTable[!train, x.var])
backward.pred<-predict(mdl.ft$backward$model,
                       newdata=baseTable[!train, x.var])
## Predicting Testset
predMat.test<-data.frame(Date=baseTable[!train, "Date"],
                    TrueValue=baseTable[!train, "PerEURO"],
                    Linear=lm.pred,
                    AICModel=aicMdl.pred,
                    BICModel=bicMdl.pred,
                    BackModel=backward.pred,
                    Ridge=ridge.pred,
                    PCR=pcr.pred,
                    PLS=pls.pred)

## Predicting Trainset
predMat.train<-data.frame(Date=baseTable[train, "Date"],
                TrueValue=baseTable[train, "PerEURO"],
                Linear=predict(mdl.ft$linear$model),
                AICModel=predict(mdl.ft$aicMdl$model),
                BICModel=predict(mdl.ft$bicMdl$model),
                BackModel=predict(mdl.ft$backward$model),
                Ridge=predict(mdl.ft$ridge$model),
                PCR=predict(mdl.ft$PCR$model, ncomp = pcr.sc),
                PLS=predict(mdl.ft$PLS$model, ncomp = pls.sc))

names(predMat.train)[grep("PCR", names(predMat.train))]<-paste("PCR.Comp", pcr.sc, sep="")
names(predMat.train)[grep("PLS", names(predMat.train))]<-paste("PLS.Comp", pls.sc, sep="")

predMat<-rbind(train=predMat.train, test=predMat.test)
stkPredMat<-melt(list(train=predMat.train, test=predMat.test), 1:2)
stkPredMat$L1<-factor(stkPredMat$L1, levels = c("train", "test"))

predMat.rpSumry<-ddply(stkPredMat, .(variable, L1), summarize,
      RMSEP=sqrt(1/length(value)*sum((TrueValue-value)^2)),
      R2pred=1-(sum((TrueValue-value)^2)/sum((TrueValue-mean(TrueValue))^2)))

## ----testPredErr, echo=FALSE---------------------------------------------
errMat<-lapply(3:ncol(predMat.test), function(x){rmserr(predMat.test[,2], predMat.test[,x])})
names(errMat)<-names(predMat.test)[-c(1:2)]
errStkMat<-melt(errMat)
errStkMat$L1<-factor(errStkMat$L1, levels = names(errMat))


## ----gofSumry, echo=FALSE------------------------------------------------
gofSumry<-ldply(names(mdl.ft)[-c(2:4)], function(x){
    data.frame(Model=x,
               AIC=AIC(mdl.ft[[x]][[2]]), 
               BIC=AIC(mdl.ft[[x]][[2]], 
                       k = log(nrow(mdl.ft[[x]][[3]]))),
               `R-Sq`=summary(mdl.ft[[x]][[2]])$r.squared,
               `R-Sq Adj`=summary(mdl.ft[[x]][[2]])$adj.r.squared,
               `Sigma`=summary(mdl.ft[[x]][[2]])$sigma,
               `F-value`=summary(mdl.ft[[x]][[2]])$fstat[1],
               `P-value`=signif(pf(summary(mdl.ft[[x]][[2]])$fstat[1], 
                            summary(mdl.ft[[x]][[2]])$fstat[2], 
                            summary(mdl.ft[[x]][[2]])$fstat[3], 
                            lower.tail = FALSE), 3))
})



## ----ValdSumry, echo=FALSE, results='hide'-------------------------------
ValdSumry<-rbind(predMat.rpSumry, data.frame(variable=rownames(cvStat), L1="cv", cvStat, row.names = NULL))
names(ValdSumry)<-c("Model", "Type", "RMSEP", "R2pred")
vs.cast<-dcast(melt(ValdSumry, 1:2), Model~Type+variable)[, c(1:3,6:7,4:5)]

ValdSumryTabl<-xtable(vs.cast, digits = 4)
caption(ValdSumryTabl)<-"Validation result containing RMSEP and R2pred for training set, cross-validation set and test set"
label(ValdSumryTabl)<-"tbl:valdSumry"
align(ValdSumryTabl)<-"lrrrrrrr"
tblHeader<-paste("\\hline Model & 
                 \\multicolumn{2}{c}{Training} & 
                 \\multicolumn{2}{c}{Cross Validation} & 
                 \\multicolumn{2}{c}{Test} \\\\ 
                 \\cline{2-7} &", 
                 paste(rep(c("RMSEP", "R2pred"), 3), 
                       collapse=" & "), 
                 '\\\\')

## ----ValdSumryPlotSetup, echo=FALSE--------------------------------------
vss<-ddply(ValdSumry, .(Type), summarize,
      Model.rmsep=Model[which.min(RMSEP)],
      Model.r2pred=Model[which.max(R2pred)],
      RMSEP=min(RMSEP),
      R2pred=max(R2pred))
vss1<-filter(melt(vss,1:3), variable=='RMSEP')[,-3]
vss2<-filter(melt(vss,1:3), variable=='R2pred')[,-2]
names(vss1)<-names(vss2)<-c("Type", "Model", "variable", "value")
vss<-rbind(vss1, vss2)





## ----coefMat, echo=FALSE-------------------------------------------------
coefMat<-cbind(sapply(c(1,4), function(x){coef(mdl.ft[[x]][[2]])[-1]}), 
               coef(mdl.ft$PCR$model, ncomp = pcr.min.comp), 
               coef(mdl.ft$PLS$model, ncomp=pls.min.comp))
coefMat<-data.frame(variable=rownames(coefMat), coefMat, row.names = NULL)
names(coefMat)<-c("vars","linear", "ridge", "pcr", "pls")
coefMat$vars<-factor(coefMat$vars, levels = coefMat$vars[order(coefMat$linear)])



## ----coefMatPrint, echo=FALSE, results='asis'----------------------------
xtable(do.call(cbind, lapply(1:2, function(x){
    rbind(coefMat[apply(coefMat[,-1], 2, order)[1:4,x],c(1,x+1)],
          coefMat[apply(coefMat[,-1], 2, order, decreasing=T)[1:4,x],c(1,x+1)])
})), digits = 4)
xtable(do.call(cbind, lapply(3:4, function(x){
    rbind(coefMat[apply(coefMat[,-1], 2, order)[1:4,x],c(1,x+1)],
          coefMat[apply(coefMat[,-1], 2, order, decreasing=T)[1:4,x],c(1,x+1)])
})), digits = 4)


## ----forecast, echo=FALSE, fig.cap="Prediction made on trained and test dataset using different models", fig.height=9.5, fig.width="\\textwidth"----
ggplot(stkPredMat, aes(Date, value))+
    geom_line(aes(color="red"))+
    facet_wrap(~variable+L1, 
               scale="free_x", 
               ncol = 4)+
    geom_line(aes(y=TrueValue, color="blue"), 
              shape=21)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45, hjust=0.5, vjust=0.5),
          text=element_text(size=9),
          legend.title=element_blank(),
          legend.position="top")+
    geom_text(data=predMat.rpSumry, 
              aes(label=paste("RMSEP:", round(RMSEP, 3), 
                              "\nR2pred", round(R2pred, 3))), 
              x=-Inf, y=Inf, hjust=-0.1, vjust=1.1, size=2.5)+
    scale_color_manual(values=c("red", "blue"), 
                       labels=c("Predicted", "Original"))




## ----chapter5-include, child="Include/Chapter-5.Rnw", eval=TRUE----------




## ----appendixVarUsed, child="Include/Appendix-varUsed.Rnw", eval=TRUE----


## ----dataDescData, echo=FALSE, warning=FALSE, results='hide'-------------
dataDescription<-read.xls(data.path, sheet = 2)


## ----dataDescTable, echo=FALSE, results='asis'---------------------------
dataDescription[,1]<-paste("\\texttt{", dataDescription[,1], "}", sep="")
names(dataDescription)[1:2]<-c("Code", "Description")
dataDescTab<-xtable(dataDescription[,1:2], align = "llX", caption = "Variable codes and their descriptions used in this paper")
print(dataDescTab, include.rownames = F, tabular.environment = "tabularx", width = "\\textwidth", floating=FALSE, booktabs = TRUE, add.to.row = list(pos = list(0),command = "\\hline \\endhead "), sanitize.text.function = function(x){x})




## ----pkgsUsed, child="Include/Appendix-pkgsUsed.Rnw", eval=TRUE----------


## ----pkgsUsed, echo=FALSE------------------------------------------------
pkgsDesc<-ldply(c(req.package, "graphics", "grDevices", "utils", "datasets", "methods", "base"), function(x){
    data.frame(
    `Package Name`=packageDescription(x)$Package,
    `Version`=packageDescription(x)$Version,
    `Title`=packageDescription(x)$Title)
})
citeKey<-c('car2011FJnWS','dplyr2014WHFR','gdata2014WG','ggplot22009WH','gridExtra2012AB','knitr2013XY','leaps2009LT','MASS2001WNV','mixlm2014SK','pls2013MBH','plyr2011WH','R2014Rcore','reshape22007WH','scales:2014Wickham','ridge2014CE','xtable2014DD','zoo2005ZAGG')
ckSrtd<-unlist(lapply(paste("^",pkgsDesc$Package.Name, sep=""), function(x){
    grep(x, x = citeKey, value = TRUE)
}))
ckSrtd<-c(ckSrtd,rep('R2014Rcore', 6))
citeCmd<-paste("\\cite{",ckSrtd,"}", sep="")






## ----revPlots, child="Include/Appendix-revPlots.Rnw", eval=TRUE----------









## ----appendixCodeUsed, child="Include/Appendix-PLSflowchart.Rnw", eval=TRUE----




## ----appendixCodeUsed, child="Include/Appendix-codeUsed.Rnw", eval=TRUE----




