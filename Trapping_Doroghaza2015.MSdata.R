library(lattice)
library(lmtest)
library(pscl)
library(MASS)
library(lme4)
library(plotrix)
library(emmeans)

t.df <- read.csv("Trapping_Doroghaza2015.csv")

nrow(t.df)
summary(t.df)
head(t.df)



d <- as.Date(as.character(t.df$Date), format="%m/%d/%y")
t.df$Date <- d
d <- as.character(t.df$Time)
d <- matrix(as.numeric(unlist(strsplit(d, ":"))), ncol=2, byrow=TRUE)
t.df$hours <- d[,1] + d[,2]/60
t.df$Grid <- as.factor(t.df$Grid)
t.df$Grid
d <- factor(as.character(t.df$Sex), labels=c("female", "male"))
t.df$Sex <- d
d <- factor(as.character(t.df$Capt), labels=c("new", "recapture"))
t.df$Capt <- d
d <- as.character(t.df$ID)
i <- d == "15" & t.df$Date == "2015-05-11"
d[i] <- "115"
t.df$ID <- factor(d)
head(t.df)
##Calculate some distance meassures

my.aggregate <- function(i, t.df) {
  #print(unique(t.df[i, "ID"]))
  l.i <- ifelse(length(i) == nrow(t.df), sum(i), length(i))
  if(l.i > 1) {
    d <- dist(t.df[i, c("trap_X", "trap_Y")])
    m <- diag(1, nrow=l.i-1)
    m <- cbind(m, rep(0, nrow(m)))
    m <- rbind(rep(0, ncol(m)), m)
    s.d <- as.matrix(d)[m > 0.5]
  } else {
    d <- NA
    s.d <- NA
  }
  a <- data.frame(
    ID=unique(as.character(t.df[i, "ID"])),
    Grid=unique(as.character(t.df[i, "Grid"])),
    sex=unique(as.character(t.df[i, "Sex"])),
    B.length=mean(t.df[i, "B_length"], na.rm=TRUE),
    B.width=mean(t.df[i, "B_width"], na.rm=TRUE),
    Tusk=mean(t.df[i, "Tusk"], na.rm=TRUE),
    days=as.numeric(max(t.df[i,"Date"]) - min(t.df[i, "Date"])),
    no.recaptures=l.i-1,
    tot.dist=sum(s.d),
    max.dist=max(d)
  )
  a
}

a.df <- data.frame()
for(id in unique(t.df$ID)) {
  i <- t.df$ID == id
  r <- my.aggregate(i, t.df)
  a.df <- rbind(a.df, r)
}
a.df$ID <- factor(a.df$ID)
a.df$sex <- factor(a.df$sex)
a.df

##Body size vs tusk size

sigmoid.fn <- function(x, x0, L.min, L.max, k) {
  L.min + (L.max - L.min)/(1 + exp(-k*(x - x0)))
}
l.sig <- nls(Tusk~sigmoid.fn(B.width, x0, L.min, L.max, k), a.df,
             start=list(x0=11.8, L.min=0.5, L.max=3.5, k=2))
l.lin <- nls(Tusk~a + b*B.width, a.df, start=list(a=-8, b=0.85))
anova(l.lin, l.sig, test="Chisq")
summary(l.sig)

####Sigmoid curve fits better than strait line does.

plot(Tusk~B.width, a.df, xlab="body width [mm]", ylab="tusk size [mm]",las=1)
rug(a.df$B.width[a.df$sex=="female"])
bwidth <- seq(min(a.df$B.width, na.rm=TRUE), 
              max(a.df$B.width, na.rm=TRUE), length.out=100)
c.sig <- coef(l.sig)
lines(bwidth, sigmoid.fn(bwidth, c.sig["x0"], c.sig["L.min"],
                         c.sig["L.max"], c.sig["k"]))
abline(v=c.sig["x0"], lty=2, col="grey")


####The rug represents the females’ body width values, The solid line is the fitted sigmoid curve, the vertical dashed line marks the inflection points.

#Analyses of movement


##Number of recaptures

na.df <- a.df[!is.na(a.df$B.width),]
na.df <- na.df[na.df$Grid == "1",] #### Only data from Grid 1 has been used! 


with(a.df, tapply(no.recaptures, sex, mean))

with(a.df, tapply(no.recaptures, sex, sd))

gof <- function(model) {
  pchisq(deviance(model), df.residual(model), lower.tail=FALSE)
}


### Scale the body width within sexes:

fem.df<-na.df[na.df$sex=="female",]
mal.df<-na.df[na.df$sex=="male",]
mean(fem.df$B.width)
sd(fem.df$B.width)
mean(mal.df$B.width)
sd(mal.df$B.width)

fem.df$rel.bodywidth <- as.numeric(scale(fem.df$B.width))
mal.df$rel.bodywidth <- as.numeric(scale(mal.df$B.width))

rel.na.df<-rbind(fem.df,mal.df)
nrow(rel.na.df)
nrow(na.df)

### Now check the models

rel.p.rec <- glm(no.recaptures ~ rel.bodywidth*sex, rel.na.df, family=poisson)
anova(rel.p.rec, test="Chisq")

rel.p.rec1 <- update(rel.p.rec, .~. - rel.bodywidth:sex)

anova(rel.p.rec, rel.p.rec1, test="Chisq")

anova(rel.p.rec1, test="Chisq")

plot(rel.p.rec1)
summary(rel.p.rec1)

gof(rel.p.rec1)

rel.qp.rec <- glm(no.recaptures ~ rel.bodywidth * sex, rel.na.df, family=quasipoisson)
anova(rel.qp.rec, test="F")

rel.qp.rec1 <- update(rel.qp.rec, .~. - rel.bodywidth:sex)
anova(rel.qp.rec, rel.qp.rec1, test="F")

anova(rel.qp.rec1, test="F")

plot(rel.qp.rec1)

summary(rel.qp.rec1)

rel.nb.rec <- glm.nb(no.recaptures ~ rel.bodywidth * sex, rel.na.df)
anova(rel.nb.rec)
drop1 (rel.nb.rec, test="Chisq")

rel.nb.rec1 <- update(rel.nb.rec, .~. - rel.bodywidth:sex)
anova(rel.nb.rec, rel.nb.rec1)

anova(rel.nb.rec1)

plot(rel.nb.rec1)

summary(rel.nb.rec1)

gof(rel.nb.rec1)

rel.h.rec <- hurdle(no.recaptures ~ rel.bodywidth * sex, rel.na.df)
summary(rel.h.rec)

rel.h.rec1 <- update(rel.h.rec, .~. - rel.bodywidth:sex)
lrtest(rel.h.rec, rel.h.rec1)

rel.h.rec2 <- update(rel.h.rec1, .~. - rel.bodywidth)
lrtest(rel.h.rec1, rel.h.rec2)

rel.h.rec2 <- update(rel.h.rec1, .~. - sex)
lrtest(rel.h.rec1, rel.h.rec2)

summary(rel.h.rec2)

#### We now compare the count models.

logLik(rel.p.rec1)

logLik(rel.nb.rec1)

logLik(rel.h.rec1)

AIC(rel.p.rec1)

AIC(rel.nb.rec1)

AIC(rel.h.rec1)

#### seems that NB fits the best

summary(rel.nb.rec1)
drop1(rel.nb.rec1,test="Chisq")

rel.nb.rec2<-update(rel.nb.rec1, .~. -rel.bodywidth)
summary(rel.nb.rec2)
drop1(rel.nb.rec2,test="Chisq")

library(plotrix)

fem.se<-std.error(rel.na.df$no.recaptures[na.df$sex=="female"], na.rm=TRUE)
fem.se

m.se<-std.error(rel.na.df$no.recaptures[na.df$sex=="male"], na.rm=TRUE)
m.se

t.se<-std.error(rel.na.df$no.recaptures, na.rm=TRUE)
t.se
recapplot <- boxplot(rel.na.df$no.recaptures~rel.na.df$sex, ylim=c(0, max(rel.na.df$no.recaptures+t.se)), xlab="Sex", 
                     ylab="Number of recaptures", cex.axis=1.1, cex=0.9, cex.lab=1.3, las=1, lwd=2,font.lab=2, font.main=4,col=c("white","black"),medcol=c("black","white") )


##Days seen

rel.l.d <- glm(days~no.recaptures*sex*rel.bodywidth, rel.na.df, family=poisson,
               subset=no.recaptures > 0)
anova(rel.l.d, test="Chisq")

summary(rel.l.d)

gof(rel.l.d)

rel.nb.days <- glm.nb(days ~ no.recaptures * sex * rel.bodywidth, rel.na.df,
                      subset=no.recaptures > 0)
anova (rel.nb.days, test="Chisq")
drop1 (rel.nb.days, test="Chisq")
gof(rel.nb.days)
AIC(rel.l.d)
AIC(rel.nb.days)

rel.nb.d1 <- update(rel.nb.days, .~. - no.recaptures:sex:rel.bodywidth)
drop1(rel.nb.d1, test="Chisq")


rel.nb.d2 <- update(rel.nb.d1, .~. - no.recaptures:sex)
drop1(rel.nb.d2, test="Chisq")

rel.nb.d3 <- update(rel.nb.d2, .~. - sex:rel.bodywidth)
drop1(rel.nb.d3, test="Chisq")

rel.nb.d4 <- update(rel.nb.d3, .~. - no.recaptures:rel.bodywidth)
drop1(rel.nb.d4, test="Chisq")


rel.nb.d5 <- update(rel.nb.d4, .~. - no.recaptures)
drop1(rel.nb.d5, test="Chisq")

rel.nb.d6 <- update(rel.nb.d5, .~. - rel.bodywidth)
drop1(rel.nb.d6, test="Chisq")

anova(rel.nb.d6)

gof(rel.nb.d6)

plot(rel.nb.d6)

summary(rel.nb.d6)

rel.mna.df<-rel.na.df[rel.na.df$no.recaptures>0,]
rel.mna.df



daysplot <- boxplot(rel.mna.df$days~rel.mna.df$sex, ylim=c(0, max(rel.mna.df$days+1)), xlab="Sex", 
                    ylab="Number of days", cex.axis=1.1, cex=0.9, cex.lab=1.3, las=1, font.lab=2,lwd=2,col=c("white","black"),medcol=c("black","white") )

## Total distance travelled


rel.mna.df<-rel.na.df[rel.na.df$no.recaptures>0,]

rel.l.td <- lm(log(tot.dist+1)~no.recaptures*sex*rel.bodywidth, rel.mna.df)
drop1(rel.l.td, test="F")

rel.l.td1 <- update(rel.l.td, .~. - no.recaptures:sex:rel.bodywidth)
drop1(rel.l.td1, test="F")

rel.l.td2 <- update(rel.l.td1, .~. - sex:rel.bodywidth)
drop1(rel.l.td2, test="F")

rel.l.td3 <- update(rel.l.td2, .~. - no.recaptures:rel.bodywidth)
drop1(rel.l.td3, test="F")
summary(rel.l.td3)

rel.l.td4 <- update(rel.l.td3, .~. - no.recaptures:sex)
drop1(rel.l.td4, test="F")
summary(rel.l.td4)

#### Plot

summary(rel.mna.df)

n.data <- data.frame(no.recaptures = 
                       rep(seq(min(rel.mna.df$no.recaptures),
                               max(rel.mna.df$no.recaptures), length=100), 2),
                     sex = rep(c("female", "male"), c(100, 100)),
                     rel.bodywidth=mean(rel.mna.df$rel.bodywidth))
pred.td <- predict(rel.l.td4, newdata=n.data, type="response", se.fit=TRUE)
n.data$tot.dist <-(pred.td$fit)
#rel.na.df$tn.se<-std.error(rel.na.df$tot.dist, na.rm=TRUE)
plot(tot.dist~no.recaptures, n.data, type="n", xlab="Number of recaptures", 
     ylab="Total distance travelled [m]",las=1, cex.axis=1.1, cex=0.9, cex.lab=1.3, las=1, font.lab=2, font.main=4, ylim=c(0, 20),xlim=c(0,10))
points(tot.dist~no.recaptures, rel.na.df, subset=sex=="female", pch=1,
       col=1)
points(tot.dist~no.recaptures, rel.mna.df, subset=c(sex=="male"), pch=17,
       col=1)

lines(exp(tot.dist) ~ no.recaptures, n.data, subset=sex=="female", lty=3,lwd=2)
lines(exp(tot.dist) ~ no.recaptures, n.data, subset=c(sex=="male" & no.recaptures<=3), col=1,lwd=2)
legend(x=0,y=21, legend=c("female", "male"),
       pch=c(1,17), lty=c(3,1), bty="n", cex=1.2,y.intersp =0.6,x.intersp = 0.1,ncol=1,lwd=2)


### Total distance travelled a day

rel.l.tdd <- lm(log((tot.dist+1)/(days+1))~sex*rel.bodywidth*no.recaptures, rel.mna.df)
drop1(rel.l.tdd, test="F")

rel.l.tdd1 <- update(rel.l.tdd, .~. - sex:rel.bodywidth:no.recaptures)
drop1(rel.l.tdd1, test="F")

rel.l.tdd2 <- update(rel.l.tdd1, .~. -sex:no.recaptures)
drop1(rel.l.tdd2, test="F")

rel.l.tdd3<-update(rel.l.tdd2, .~. -sex:rel.bodywidth)
drop1(rel.l.tdd3, test="F")

rel.l.tdd4<-update(rel.l.tdd3, .~. -rel.bodywidth:no.recaptures)
drop1(rel.l.tdd4, test="F")
summary(rel.l.tdd4)

rel.l.tdd5<-update(rel.l.tdd4, .~. -rel.bodywidth)
drop1(rel.l.tdd5, test="F")

rel.l.tdd6<-update(rel.l.tdd5, .~. -no.recaptures)
drop1(rel.l.tdd6, test="F")

rel.mna.df$dayly<-(rel.mna.df$tot.dist+1)/(rel.mna.df$days+1)
summary(rel.mna.df$dayly)
daysplot <- boxplot(rel.mna.df$dayly~rel.mna.df$sex, ylim=c(0, 5), xlab="Sex", 
                    ylab="Distance travelled per day [m]", cex.axis=1.1, cex=0.9, cex.lab=1.3, las=1, font.lab=2,lwd=2,col=c("white","black"),medcol=c("black","white"))

#### Plot for short comm

par(mfrow = c(1, 2))

recapplot <- boxplot(rel.na.df$no.recaptures~rel.na.df$sex, ylim=c(0, max(rel.na.df$no.recaptures+t.se)), xlab="Sex", main="(a)",
                     ylab="Number of recaptures", cex.axis=1.1, cex=0.9, cex.lab=1.3, las=1, lwd=2,font.lab=2, font.main=4,col=c("white","black"),medcol=c("black","white") )

daysplot <- boxplot(rel.mna.df$days~rel.mna.df$sex, ylim=c(0, max(rel.mna.df$days+1)), xlab="Sex", main="(b)",
                    ylab="Number of days", cex.axis=1.1, cex=0.9, cex.lab=1.3, las=1, font.lab=2,lwd=2,col=c("white","black"),medcol=c("black","white") )


