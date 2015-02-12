
## ----frontMatter, child="frontMatter.Rnw"--------------------------------


## ----LoadingPkgs, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
req.package<-c("pls", "xtable", "MASS", "car", "corrplot", "gdata", "dplyr", "ggplot2", "reshape2","ridge", "grid", "gridExtra", "stargazer", "devtools", "ggbiplot", "mixlm", "Hmisc", "knitr")
lapply(req.package, require, character.only=TRUE, quietly = T, warn.conflicts = F)




## ----setup, include=FALSE, cache=FALSE, echo=FALSE-----------------------
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
plotTS<-function(dataSet, vars, colVar=NA){
  plt<-ggplot(filter(dataSet, variable%in%vars), aes(Date))
  ifelse(is.na(colVar),
  plt<-plt+geom_line(aes(y=value), size=.3),
  plt<-plt+geom_line(aes_string(y="value", col=colVar), size=.3)    
  )
  plt<-plt+theme_bw()
  plt<-plt+theme(axis.text=element_text(size=7), strip.text.x=element_text(size=8))
  plt<-plt+facet_wrap(~variable, scales = "free_y", ncol=3)
  plt<-plt+theme(legend.position="top", legend.title=element_blank())
  return(plt)
}

## Plotting Model Coefficients with their state of significance
test.plot<-function(coef.matrix, alpha=0.05){
  .e<-environment()
  names(coef.matrix)<-c("Estimate", "StdError", "t.value", "p.value")
  idx<-order(row.names(coef.matrix))
  cp<-ggplot(coef.matrix[idx,], aes(x=row.names(coef.matrix[idx,]), y=t.value), environment = .e)
  cp<-cp+geom_bar(stat="identity", position = "identity",
                  fill=ifelse(coef.matrix[idx,"p.value"]<alpha, "firebrick3", "dodgerblue3"))
  cp<-cp+geom_text(aes(y=ifelse(coef.matrix[idx, "t.value"]>0,t.value+0.5, t.value-0.5), 
                       label=round(coef.matrix[idx,"Estimate"], 2)), angle=90, size=3)
  cp<-cp+theme_bw()+labs(x="", y="T-Value")
  cp<-cp+theme(axis.text.x=element_text(angle=90, hjust=1, size=10))
  cp<-cp+scale_fill_manual("Status", values=c("firebrick3", "dodgerblue3"), 
                           labels=c("Significant", "Non-Significant"))
  cp<-cp+geom_hline(yintercept=c(-1,1)*qt(alpha/2, df = abs(diff(dim(coef.matrix))), lower.tail = F), 
                    color="red", linetype="dashed")
  cp<-cp+theme(legend.title=element_blank(), 
               legend.position=c(0.8,0.9))
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
diagPlot<-function(model){
  p1<-ggplot(model, aes(.fitted, .resid))+geom_point()
  p1<-p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
  p1<-p1+xlab("Fitted values")+ylab("Residuals")
  p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_bw()
  
  p2<-ggplot(model, aes(qqnorm(as.vector(.stdresid))[[1]], .stdresid))+geom_point(na.rm = TRUE)
  p2<-p2+geom_abline(aes(qqline(as.vector(.stdresid))))+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
  p2<-p2+ggtitle("Normal Q-Q")+theme_bw()
  
  p3<-ggplot(model, aes(.fitted, sqrt(abs(.stdresid))))+geom_point(na.rm=TRUE)
  p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
  p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
  p3<-p3+ggtitle("Scale-Location")+theme_bw()
  
  p4<-ggplot(model, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity")
  p4<-p4+xlab("Obs. Number")+ylab("Cook's distance")
  p4<-p4+ggtitle("Cook's distance")+theme_bw()
  
  p5<-ggplot(model, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm=TRUE)
  p5<-p5+stat_smooth(method="loess", na.rm=TRUE)
  p5<-p5+xlab("Leverage")+ylab("Standardized Residuals")
  p5<-p5+ggtitle("Residual vs Leverage Plot")
  p5<-p5+scale_size_continuous("Cook's Distance", range=c(1,5))
  p5<-p5+theme_bw()+theme(legend.position="bottom")
  
  p6<-ggplot(model, aes(.hat, .cooksd))+geom_point(na.rm=TRUE)+stat_smooth(method="loess", na.rm=TRUE)
  p6<-p6+xlab("Leverage hii")+ylab("Cook's Distance")
  p6<-p6+ggtitle("Cook's dist vs Leverage hii/(1-hii)")
  p6<-p6+geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")
  p6<-p6+theme_bw()
  
  return(list(rvfPlot=p1, qqPlot=p2, sclLocPlot=p3, cdPlot=p4, rvlevPlot=p5, cvlPlot=p6))
}

## Generate summary plot from a fitted model to annotate other plot
sumryBlock<-function(model){
  return(paste("R-Sq = ",signif(summary(model)$r.squared, 2),
                       "\nAdj R-Sq =",signif(summary(model)$adj.r.squared, 2),
                       "\nSigma =",signif(summary(model)$sigma, 2),
                       "\nF =",signif(as.vector(summary(model)$fstatistic[1]), 2))
         )
}

model.sumry<-function(model, call=TRUE, coefMat=TRUE, sumry=TRUE){
    if("lm"%nin%class(model)){
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
            cat("\n sigma: ",s, " on ", df[2], " degrees of freedom\n",
            "Multiple R-squared: ", r.sq, "\n",
            "Adjusted R-Squared: ", adj.r.sq,"\n",
            "F-statistic: ", f, "on", f.df.num,"and",f.df.den, "degree of freedom\n",
            "P-value: ", pf(f, f.df.num, f.df.den, lower.tail = FALSE))
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




## ----data-prep, child="Include/DataPreperation.Rnw"----------------------


## ----dataSetup, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
baseTable<-read.xls(data.path, sheet = "FinalData")
baseTable[,1]<-as.Date(baseTable[,1], format="%d/%m/%Y")
baseTable[,"Testrain"]<-as.logical(baseTable[,"Testrain"])


## Label Variables in baseTable
labelTable<-read.xls(data.path, sheet = "FinalCodeBook", stringsAsFactors=FALSE)
for(i in 1:ncol(baseTable)){
    Hmisc::label(baseTable[,i])<-labelTable[i,2]
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




## ----chapter1-include, child="Include/Chapter-1.Rnw", eval=FALSE---------
## NA


## ----chapter2-include, child="Include/Chapter-2.Rnw", eval=FALSE---------
## NA


## ----chapter3-include, child="Include/Chapter-3.Rnw", eval=TRUE----------


## ----mdlFitCriteriaPlot, child="mdlFitCriteriaPlot.Rnw"------------------






## ----chapter4-include, child="Include/Chapter-4.Rnw", eval=FALSE---------
## NA


## ----commons, child="Include/Commons.Rnw"--------------------------------


## ----modelFitting, echo=FALSE, results='hide'----------------------------
pls.options(plsralg="oscorespls")
mdl<-c("lm", "pcr", "plsr", "linearRidge")
mdl.ft<-list()

for(i in seq_along(mdl)){
    mdl.ft[[i]]<-fit.model(mdl[i], y.var, x.var, baseTable, scaling=mdl[i] %in% c("plsr", "pcr"))
}

names(mdl.ft)<-c("linear", "PCR", "PLS", "ridge")

## Forward Selection Model (criteria: level of significance)
mdl.ft$lm.forward.alpha<-list(formula=mdl.ft$linear$formula, model=forward(lm(formula = mdl.ft$linear$formula, data=mdl.ft$linear$data), alpha = 0.05, full = FALSE), data=mdl.ft$linear$data)

## Backward Elimination Model (criteria: level of significance)
mdl.ft$lm.back.alpha<-list(formula=mdl.ft$linear$formula, model=backward(lm(formula = mdl.ft$linear$formula, data=mdl.ft$linear$data), alpha = 0.05, full = FALSE, hierarchy = TRUE), data=mdl.ft$linear$data)

## Forward Selection Model (criteria: AIC)
mdl.ft$lm.forward.aic<-list(formula=mdl.ft$linear$formula, model=step(lm(formula = mdl.ft$linear$formula, data=mdl.ft$linear$data), direction = "forward", trace=FALSE, k = 2), data=mdl.ft$linear$data)

## Backward Elimination Model (criteria: AIC)
mdl.ft$lm.back.aic<-list(formula=mdl.ft$linear$formula, model=step(lm(formula = mdl.ft$linear$formula, data=mdl.ft$linear$data), direction = "backward", trace=FALSE, k = 2), data=mdl.ft$linear$data)

mdl.labels<-c("Linear Model", "Principal Component Regression", "Partial Least Square Regression", "Ridge Regression", "Forward Selection Model (criteria:level of significance)", "Backward Elimination Model (criteria: level of significance)", "Forward Selection Model (criteria:AIC)", "Backward Elimination Model (criteria: AIC)")
for(i in 1:length(mdl.ft)){
    Hmisc::label(mdl.ft[[i]][[2]])<-mdl.labels[i]
    class(mdl.ft[[i]][[2]])<-rev(class(mdl.ft[[i]][[2]]))
}



## ----chapter4a-include, child="Include/Chapter-4a.Rnw", eval=FALSE-------
## NA


## ----chapter4b-include, child="Include/Chapter-4b.Rnw", eval=FALSE-------
## NA


## ----chapter4c-include, child="Include/Chapter-4c.Rnw", eval=FALSE-------
## NA


## ----appendixA, child="Include/AppendixA.Rnw"----------------------------


## ----dataDescData, echo=FALSE, warning=FALSE, results='hide'-------------
dataDescription<-read.xls(data.path, sheet = 2)


## ----dataDescTable, echo=FALSE, results='asis'---------------------------
dataDescription[,1]<-paste("\\texttt{", dataDescription[,1], "}", sep="")
names(dataDescription)[1:2]<-c("Code", "Description")
dataDescTab<-xtable(dataDescription[,1:2], align = "llX", caption = "Variable codes and their descriptions used in this paper")
print(dataDescTab, include.rownames = F, tabular.environment = "tabularx", width = "\\textwidth", floating=FALSE, booktabs = TRUE, add.to.row = list(pos = list(0),command = "\\hline \\endhead "), sanitize.text.function = function(x){x})




