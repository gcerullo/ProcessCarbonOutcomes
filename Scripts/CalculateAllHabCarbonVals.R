#09.06.24
#This code calculates carbon (above and belowground) of each habitat type through time 

#Code notes:
#Plantation ACD is calculate is from SSB inventory data 
#Primary, 1L and R is from Philipson and uses edited code originally written by 
#Philipson 2020; Science. (I adjust Philipson's code to allow estimates for 60 years (they
#estimate out to 30/35 years; I assume that slope and intercept of their models stay the same 
#and that 1L and R values plateu once they reach primary forest values)
#I estimate ACD in twice-logged forest, and use three 2L typologies 

#This code also
#1. Adds different harvest delays 
#2 Incorporates belowground carbon and necromass processes

#Code output
##NB this code outputs a master csv called "allHabCarbon_60yrACD_withDelays.csv"
# where ACD refers only to the raw above ground carbon density values and 
#where full_carbon incorporate belowground processes.


library(lme4)
library(lattice)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(broom)
library(tidyverse)
library(data.table)

#DEFINE KEY PARAMS ####
#define amount of ACD loss from second rotation (ACD loss from getting 31.2m3 of timber in second harvest)



#RATIONALE - BEING USED for 2L -> 2L: assume forest that starts out as twice logged was logged fiften years 
#after first logging (e.g. first harvest = t= -15, second harvest t = 0; tracks reality)
#1. For Chris Philipson's model of ACD recovery after once-logging,  it looks like the first harvesting rotation led to a decline of 155.9869 ACD_ha-1  (1 yrs after once-logged = 44.01312 +/- 31.20968), down  from ~200 ACD_ha-1 for primary forest, and removing 112.96 m3 (+/- 22.42) of timber in first rotation  
#2 Thus 155.9869/112.96 =  1.380904 decline in ACD per m3 harvest
#3. Again, from Philipson's models, once-logged forest has an ACD of 84.56991 (+/- 30.47371) after 15 years of recovery.  
#4. Twice-logging removes 31.24 m3 of timber (+/- 10.4) 
#5. I assume the same ACD loss per m3 harvested as once-logged 
#6. Thus 31.24 * 1.380904 = 43.13944 = ACD loss from second harvest
#7. Thus given a 15yr once-logged ACD before logging, then:
#      84.56991 - 43.94059 = 40.62932 in twice-logged in yr 0   

harvest2ndACD_loss = 43.13944
#define ACD left after second rotation (Mg.ha-1), 30 years after 1st logging 
ACD_2L_yr0_30yrAfter1L_30yrAfter1L <-84.88416  # 30y3 1L =  128.02360   -  harvest2ndACD_loss
#define ACD after for scenarios beginning as 2L
ACD_2L_starting2L_15yrAfter1L <- 40.62932

#-----read in inputs showing habitat transitions --------
hab_by_year <- read.csv("Inputs/HabByYears.csv", strip.white = TRUE) %>%  
  rename(true_year = year, 
         functionalhabAge = functional_habAge, 
         habitat = transition_habitat) %>% select(-X) %>%  
  #remove all improved plantation yielding varieties (can add in if interested to)
  filter(!str_detect(habitat, "improved"))

  



#--------- read in plantation data ---------

plantation <- read.csv("Outputs/plantation_carbon.csv") %>% select(-X)
names(plantation)


#----- carbon in primary, restored and once-logged --------
# get carbon of restored, primary and once-logged from Philipson 2020
#uses Philipson 2020 code from: https://github.com/PhilipsonChristopher/CarbonRecovery/blob/master/Fig1/Figure1Code.R
data <- read.csv("RawData/Philipson20_PlotData2.csv")

## Seperate Logged forest plots from unlogged plots
Logged <- subset(data, Forest=="Logged")
UnLogged <- subset(data, Forest=="UnLogged")

length(levels(Logged$Plot))


head(Logged)
summary(Logged)
str(Logged)

m1 <- lmer(ACD~  YearsSinceLogging  * FACE + (1|Plot)  + (1|LoggingMethod:Coupe) , REML=1 , data=Logged , na.action=na.fail) 
dotplot(ranef(m1, condVar=T))$Plot
dotplot(ranef(m1, condVar=T))$'LoggingMethod:Coupe'
anova(m1)
summary(m1)
confint(m1)

# Bootstrap errors 
## CI on restored slope
Logged$FACE <- factor(Logged$FACE)
Logged$FACE <- relevel(Logged$FACE, ref="ProjectScenario")
levels(Logged$FACE)
m1 <- lmer(ACD~  YearsSinceLogging  * FACE + (1|Plot)  + (1|LoggingMethod:Coupe) , REML=1 , data=Logged , na.action=na.fail) 
round(fixef(m1),1) 
m1_mer <- bootMer(m1, nsim=1000, FUN=fixef, ncpus=8)
round(apply(m1_mer$t, 2, quantile, c(0.025, 0.975)), 2)

## CI on Natural regen slope
Logged$FACE <- relevel(Logged$FACE, ref="Baseline")
factor(Logged$FACE)
m1 <- lmer(ACD~  YearsSinceLogging  * FACE + (1|Plot)  + (1|LoggingMethod:Coupe) , REML=1 , data=Logged , na.action=na.fail) 
round(fixef(m1),1) 
m1_mer <- bootMer(m1, nsim=1000, FUN=fixef, ncpus=8)
round(apply(m1_mer$t, 2, quantile, c(0.025, 0.975)), 2)


#  Variance attributable to plot size
#  adding a random effect for each dataset accounts for variance attributable to plot size
#  adding '+(1| Dataset)' to FullModel above does not alter the conclusions from the model (recovery slope is the same)
#  m1_b <- lmer(ACD ~  YearsSinceLogging  * FACE + (1|Plot)  + (1|LoggingMethod:Coupe) +(1| Dataset), data=Logged , na.action=na.fail) 
#  summary(m1_b)
#  summary(m1)
#  dotplot(ranef(m1_b, condVar=T))$Dataset



### 1. Plot predictions from model (m1)

# (NEWDATa<-expand.grid(YearsSinceLogging =seq(from=3, to=30, by=1), FACE=levels(Logged$FACE)[1]))
# (NEWDATb<-expand.grid(YearsSinceLogging =seq(from=8, to=35, by=1), FACE=levels(Logged$FACE)[2]))
# NEWDAT <- rbind(NEWDATa, NEWDATb)
# NEWDAT$ACD<-0

# ------Modify philipson code to etend to 60 yr estimation --------
#GC CODE; EXTEND above code TO PREDICT 60 YEARS INTO THE FUTURE 
(NEWDATa<-expand.grid(YearsSinceLogging =seq(from=0, to=60, by=1), FACE=levels(Logged$FACE)[1]))
(NEWDATb<-expand.grid(YearsSinceLogging =seq(from=0, to=60, by=1), FACE=levels(Logged$FACE)[2]))
NEWDAT <- rbind(NEWDATa, NEWDATb)
NEWDAT$ACD<-0


mm <- model.matrix(terms(m1),NEWDAT)
length(fixef(m1))
NEWDAT$ACD <- mm %*% fixef(m1)
pvar1 <- diag(mm %*% tcrossprod(vcov(m1),mm))
(NEWDAT <- data.frame(
  NEWDAT
  , plo = NEWDAT$ACD-2*sqrt(pvar1)
  , phi = NEWDAT$ACD+2*sqrt(pvar1)))
head(NEWDAT)


NEWDAT$FACE <- factor(NEWDAT$FACE)
levels(NEWDAT$FACE)[levels(NEWDAT$FACE)=="Baseline"] <- "A: Natural regeneration"
levels(NEWDAT$FACE)[levels(NEWDAT$FACE)=="ProjectScenario"] <- "B: With active restoration"

levels(Logged$FACE)[levels(Logged$FACE)=="Baseline"] <- "A: Natural regeneration"
levels(Logged$FACE)[levels(Logged$FACE)=="ProjectScenario"] <- "B: With active restoration"

levels(NEWDAT$FACE)
levels(Logged$FACE)


## add line indicating timing of restoration
RestorLine <- data.frame(x=2, y=399, xend=23, yend=399, FACE="B: With active restoration")

Fig1 <- ggplot(data= Logged, aes(YearsSinceLogging, ACD, FACE))+
  xlab("Years Since Logging")+
  ylab(expression(paste("Aboveground Carbon Density", " ","(Mg h", a^-1, sep = "", ")")))+
  scale_y_continuous(limits=c(0,399))+
  facet_grid(~FACE)+ # , scales="free_x", space="free"
  # geom_point(col="grey25", pch="O")+ # grey30 grey20
  geom_line(data= Logged,aes(YearsSinceLogging,ACD, group=Plot),  size=0.5, colour= "grey")+ # lines per plot
  geom_line(data=NEWDAT,aes(YearsSinceLogging, ACD, group=FACE),  linetype=1, size=1.2, colour= "blue")+
  geom_line(data=NEWDAT,aes(YearsSinceLogging,plo, group=FACE), size=1, linetype =3, colour= "blue")+
  geom_line(data=NEWDAT,aes(YearsSinceLogging,phi, group=FACE),size=1, linetype =3,  colour= "blue")+
  theme_bw()+ theme(strip.text.y = element_text( size = 12, hjust = .5), strip.text.x = element_text( size = 12, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(strip.text.x = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(strip.text.y = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(axis.title.y = element_text(size = 14, angle = 90))+
  theme(axis.title.x = element_text(size = 14, angle = 00))+
  theme(axis.text = element_text(size = 14, colour = "black"))+labs(title="")+
  geom_segment(data= RestorLine, mapping=aes(x=x, y=y, xend=xend, yend=yend),inherit.aes=FALSE, arrow=arrow(angle=65, ends="both", length=unit(0.1, "inches")), size=1, color="darkgreen")
Fig1


#### Single panel version (easier to compare slopes)
Fig_S2 <- ggplot(data= Logged, aes(YearsSinceLogging, ACD, FACE))+
  xlab("Years Since Logging")+
  ylab(expression(paste("Aboveground Carbon Density", " ","(Mg h", a^-1, sep = "", ")")))+
  scale_y_continuous(limits=c(0,399))+
  geom_line(data= Logged,aes(YearsSinceLogging,ACD, group=Plot, colour =FACE),  size=0.5, alpha=0.2)+ # lines per plot
  geom_line(data=NEWDAT,aes(YearsSinceLogging, ACD, colour =FACE),  linetype=1, size=1.5)+
  geom_line(data=NEWDAT,aes(YearsSinceLogging,plo, colour =FACE), size=1, linetype =2)+
  geom_line(data=NEWDAT,aes(YearsSinceLogging,phi, colour =FACE),size=1, linetype =2)+
  theme_bw()+ theme(strip.text.y = element_text( size = 12, hjust = .5), strip.text.x = element_text( size = 12, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(strip.text.x = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(strip.text.y = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(axis.title.y = element_text(size = 14, angle = 90))+
  theme(axis.title.x = element_text(size = 14, angle = 00))+
  theme(axis.text = element_text(size = 14, colour = "black"))+labs(title="")+
  scale_color_manual(values=c("blue", "red"))+ 
  theme(legend.position = "none")
Fig_S2


## Unlogged forest for comparison
head(UnLogged)
Prime <- lmer(ACD ~ 1 + (1|MeasureTime), data= UnLogged)
Prime
Prime_mer <- bootMer(Prime, nsim=1000, FUN=fixef, ncpus=8)
summary(Prime)
round(fixef(Prime),0)
round(apply(Prime_mer$t, 2, quantile, c(0.025, 0.975)), 0)
(PrimeCIs <- round(apply(Prime_mer$t, 2, quantile, c(0.025, 0.975)), 0))


### Unlogged / Primary forest
##### right hand panel (carbon potential)
head(UnLogged)
UnLogged$FACE <- "Unlogged"
PrimaryReference <- ggplot(data= UnLogged, aes(YearsSinceLogging, ACD, FACE))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits=c(0,399), labels =NULL)+  
  geom_point(col="grey25", pch="O")+
  geom_segment(mapping=aes(x=0, y=PrimeCIs[1], xend=0, yend=PrimeCIs[2]), arrow=arrow(angle=90, ends="both", length=unit(0.1, "inches")), size=0.5, color="blue") + 
  geom_point(x=0, y=round(fixef(Prime),0),  col="blue", pch="O",size=5)+
  theme_bw()+ theme(strip.text.y = element_text( size = 12, hjust = .5), strip.text.x = element_text( size = 12, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(strip.text.x = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(strip.text.y = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
  theme(axis.title.y = element_text(size = 14, angle = 90))+
  theme(axis.title.x = element_text(colour="white", size = 14, angle = 00))+
  theme(axis.text = element_text(size = 14, colour = "white"))+
  theme(axis.ticks.x=element_line(colour = "white"))+
  theme(panel.grid.minor.x = element_blank())+
  labs(title="")+ theme(panel.grid.major.x = element_blank())

(PrimaryReferenceF1 <- PrimaryReference + facet_grid(~FACE))

# ---- ACD and CI by habitat type ------
#BUILD A DATAFRAME OF ACD and confidence intervals by age for each habitat type  

#WE WILL USE PRIMARY FROM PHILIPSON (and not 50ha Jucker plot data) AS IT'S MODELLED against TIME AND ALSO HAS BOOTSTRAP MEASURES;
#AND THEREFORE IS MORE IN-LINE WITH RESTORED AND ONCE-LOGGED VALUES THAT I AM USING 

# 

#primary forest CIs 
PrimeCIs <- as.data.frame(PrimeCIs)
#primaryforest mean values 
mean_primary_val <- round(fixef(Prime),0)

primary_Vals <- data.frame(ACD = mean_primary_val, 
                           lwr_ACD = PrimeCIs[1,], 
                           upr_ACD = PrimeCIs[2,],
                           habitat = "primary")
rownames(primary_Vals) <- NULL

#get max primary forest val
max_val <- primary_Vals$ACD %>% unique

# -------------get once logged and restored mean ACD and 95 CI by age -------------

L1_R <- NEWDAT %>% 
  rename(habitat = FACE,
         lwr_ACD = plo,
         upr_ACD = phi) %>%  
  mutate(habitat = case_when(
    habitat == "A: Natural regeneration" ~ "once-logged",
    habitat == "B: With active restoration" ~ "restored")) %>%  
  rename(time_since_intervention = YearsSinceLogging) %>%  
  mutate(true_age = time_since_intervention) %>%  
  #MANUALLY CAUSE CARBON TO PLATEU IN REGENERATING FOREST ONCE IT REACHES OG VALUES 
  mutate(
    ACD = ifelse(ACD >max_val , NA, ACD),
    lwr_ACD = ifelse(ACD > max_val, NA, lwr_ACD),
    upr_ACD = ifelse(ACD > max_val, NA, upr_ACD)) %>%  
  #take on the max value if NA
  mutate(
    ACD = ifelse(is.na(ACD), max(ACD, na.rm = TRUE), ACD),
    lwr_ACD = ifelse(is.na(lwr_ACD), max(lwr_ACD, na.rm = TRUE), lwr_ACD),
    upr_ACD = ifelse(is.na(upr_ACD), max(upr_ACD, na.rm = TRUE), upr_ACD)
  )

oncelogged <- L1_R %>% filter(habitat == "once-logged") %>% filter(true_age < 31)

#------------- get correct years -------------
#make corrections so that we have values for all (and for the correct) years


#----60 yrs of plantations ----

#1.make plantation carbon over 60 years
plantation <-  plantation %>% 
  rename(time_since_intervention = plantationAge)  

ec_plantation_grid <- data.frame(true_age = seq(1,60, by = 1), 
                                 time_since_intervention = seq(1,6, by = 1))
al_plantation_grid <- data.frame(true_age = seq(1,60, by = 1), 
                                 time_since_intervention = seq(1,12, by = 1))


#complete cases of plantation by age
plantation_ec <- plantation %>% filter(habitat == "eucalyptus") %>% 
  right_join(ec_plantation_grid, by = "time_since_intervention") 
plantation_al <- plantation %>% filter(habitat == "albizia") %>% 
  right_join(al_plantation_grid, by = "time_since_intervention") 

plantation <- rbind(plantation_ec, plantation_al)

#add in zero plantation age (for the very first year of clearance)
zero_eucPlant <- data.frame(true_age = 0,
                            ACD = 0, 
                            lwr_ACD = 0, 
                            upr_ACD = 0, 
                            habitat = "eucalyptus", 
                            time_since_intervention = 0)
zero_albPlant <- data.frame(true_age = 0,
                            ACD = 0, 
                            lwr_ACD = 0, 
                            upr_ACD = 0, 
                            habitat = "albizia", 
                            time_since_intervention = 0)


plantation <- zero_eucPlant %>% rbind(zero_albPlant) %>% rbind(plantation)


#currently have eucalyptus and albizia carbon (not for improved, or current yields, which will actually probs have diff carbon!)
plantation <- plantation %>% mutate(habitat = case_when(
  habitat == "eucalyptus" ~ "eucalyptus_current",
  habitat == "albizia" ~ "albizia_current",
  TRUE ~ habitat
))

#Uncommented if interested in adding info on improved plantation variety 
#improved_ec <- plantation %>% filter(habitat == "eucalyptus_current") %>% mutate(habitat = "eucalyptus_improved")
#improved_al <- plantation %>% filter(habitat == "albizia_current") %>% mutate(habitat = "albizia_improved")
#plantation <- plantation %>% rbind(improved_al) %>% rbind(improved_ec)

#----60 yrs of primary ----
#2.make primary carbon over 60 years

primary_Vals <- data.frame(true_age = seq(1,60, by = 1)) %>% 
  cbind(primary_Vals[rep(1, 60), ]) %>% 
  mutate(time_since_intervention = true_age)

#add zero year
zeroPrim <- data.frame(true_age = 0,
                       ACD = 203, 
                       lwr_ACD = 160, 
                       upr_ACD = 250, 
                       habitat = "primary", 
                       time_since_intervention = 0)

primary_Vals <- zeroPrim %>% rbind(primary_Vals)


#----60 yrs of twice-logged ----
#--- three 2L typologies (1,2,3) ----

#Briefly: 
#1 = P -> 2L
#2 = 2L -> 2L 
#3 = 1L -> 2L 


#NOTE THERE ARE THREE TYPES OF TWICE-LOGGING CURVES THRU TIME BUILT HERE: 
#1. ACD_2L_yr0_30yrAfter1L_30yrAfter1L - second rotations happens 30 years after 1st logging. Always the case in all_primary and all_primary_deforested scenarios

#2, ACD_2L_starting2L_15yrAfter1L - second rotation happens 15 years after 1st rotation (in t-1) (as in reality) 
#- this is the case for all scenarios that start 2L

#and also is true for scenarios that go from 1L (at the beginning of the scenarios) -> 2L....
#....which is operationalised by not allowing 1L-2L transitions for 15 years....see below

#3. twice_logged_delays  - second rotation happens at varying intervals after the once-logging from 15-30 years, depending on the harvesting delay. 
# this is important for scenarios that begin with once-logged forest (eg. mostly_1L and mostly_1L_deforested scenarios) 
# the delay matters because if the second harvest happens after a big time delay (e.g. 25 yrs) then ACD was recovering all this time in once-logged, so the 2nd harvest leaves behind higher ACD 
#NB - no logging of 1L -> 2L is allowed in the first 15 yrs, as the forest was logged in t-1. 


#----#calculate for 1&2 (twice-logged) ------

#Above ground carbon in twice-logged straight after logging (y-intercept)
print(ACD_2L_yr0_30yrAfter1L_30yrAfter1L)  #this assumes second logging happens 30 years after once-logging 
print(ACD_2L_starting2L_15yrAfter1L) #this assume second logging happens 15 yrs after first logging, at year t-1...this is the ACD for scenarios starting as twice logged

# Define the y-intercept, slope, and confidence intervals
Slope_1L_df <- L1_R %>% filter(habitat == "once-logged")   
Slope_1L <-  lm(ACD ~ time_since_intervention, data = Slope_1L_df)
# Extract the slope using tidy()
Slope_1L <- tidy(Slope_1L)$estimate[2]  
years <- 0:60

# Calculate the predicted y-values for each year
predicted2L_values <- ACD_2L_yr0_30yrAfter1L_30yrAfter1L + Slope_1L * years # #this assumes second logging happens 30 years after once-logging
predicted2L_values_starting <- ACD_2L_starting2L_15yrAfter1L + Slope_1L * years#this assume second logging happens 15 yrs after first logging, at year t-1


# Calculate the lower and upper limits of the confidence intervals
error_2L <- L1_R %>% filter(habitat == "once-logged") %>%
  mutate(error95 = ACD - lwr_ACD) %>%
  select(error95)

twice_log_fun <- function(predVals){
  
  predVals %>%  data.frame() %>%
    rename(ACD = 1) %>%  
    cbind(habitat = "twice-logged") %>%
    cbind(years = 0:60) %>% 
    cbind(error_2L) %>% 
    mutate(lwr_ACD =  ACD - error95, 
           upr_ACD =  ACD + error95) %>%  
    
    #MANUALLY CAUSE CARBON TO PLATEU IN  FOREST ONCE IT REACHES OG VALUES 
    mutate(
      ACD = ifelse(ACD >max_val , NA, ACD),
      lwr_ACD = ifelse(ACD > max_val, NA, lwr_ACD),
      upr_ACD = ifelse(ACD > max_val, NA, upr_ACD)) %>%  
    #take on the max value if NA
    mutate(
      ACD = ifelse(is.na(ACD), max(ACD, na.rm = TRUE), ACD),
      lwr_ACD = ifelse(is.na(lwr_ACD), max(lwr_ACD, na.rm = TRUE), lwr_ACD),
      upr_ACD = ifelse(is.na(upr_ACD), max(upr_ACD, na.rm = TRUE), upr_ACD)
    ) %>% select(-error95) %>%  
    rename(true_age = years) %>%  
    mutate(time_since_intervention = true_age)
  
}


error_2L <- error_2L %>% cbind(true_year = 0:60)



# for 1: P -> 2L
twice_L_30yrPost <-  twice_log_fun(predicted2L_values) %>% cbind(original_habitat = "primary")%>% select(-true_age)

#adjust primary -> 2L so that is undergoes a once-logged conversion, then a twice-logged conversion after 30 yrs 
names(twice_L_30yrPost)

yr0_primary_2L <- primary_Vals  %>% slice(1) %>% rename(true_year = true_age)  %>% 
  mutate(habitat = "twice-logged",
         original_habitat = "primary", 
         functional_habitat = "primary" )

yr1_30_primary_2L <- oncelogged  %>% rename(true_year = true_age) %>%
  filter(true_year >= 1, true_year <= 29) %>%  
  mutate(habitat = "twice-logged",
         original_habitat = "primary", 
         functional_habitat = "once-logged")

beginning_habs_primary_2L <- rbind(yr0_primary_2L,yr1_30_primary_2L) %>% 
  rename(functionalhabAge = time_since_intervention)

twice_L_30yrPost <- twice_L_30yrPost %>% mutate(true_year = true_year + 30) %>%
  rename(functionalhabAge = time_since_intervention) %>% 
  filter(true_year <61) %>% 
  mutate(functional_habitat = "twice-logged")

twice_L_30yrPost <- rbind(twice_L_30yrPost,beginning_habs_primary_2L) 

#for 2L -> 2L 
twice_L_15yrPost <-  twice_log_fun(predicted2L_values_starting) %>% cbind(original_habitat = "twice-logged") %>%  
  cbind(functional_habitat = "twice-logged") %>% 
  rename(functionalhabAge = time_since_intervention) %>% 
  select(-true_age)


from2L_2L_AND_Primary_2L <- rbind(twice_L_30yrPost) %>%  #"Recovering where yr 0 harvests 30-yr 1L"),
  rbind(twice_L_15yrPost)  # "Recovering where yr 0 harvests once AND twice")

#-------------  #calculate for 3  ----------------

oncelogged <- L1_R %>% filter(habitat == "once-logged") %>% filter(true_age <30)

#calcaulate carbon for parcel, starting with a once-logged forest, which changes 
#depending on the delay before re-harvest (as ACD is recovering in 1L forest pre-harvest)
twice_logged_delays <- oncelogged %>% mutate(harvest2ndACD_loss = harvest2ndACD_loss, 
                                             harvest_delay = time_since_intervention) %>% 
  #calculate ACD in year 0 from 2nd harvest, which is contingent on the delay time 
  mutate(ACD_twicelogged_yr0 = ACD - harvest2ndACD_loss) %>%  
  filter(harvest_delay >14) 


#only estimate 1L -> 2L after a delay of 15 years, because not allowed to harvest once-logged that harvested in t-1 imeddiately, have to wait for 15 yrs
years <-0:60
harvest_delay_vt = 15:29
rep_delay <- length(years)

# Define a function to calculate over 60 years the ACD of twice logging 1L forests after different delays 
twice_logged_delayed_fun <- function(i) {
  result <- i + Slope_1L * years
  result_df <- data.frame(ACD = result, true_year = years)
  return(result_df)
}

#estimate ACD of going from 1L -> with different time delays..
predicted2L_values_list <- lapply(twice_logged_delays$ACD_twicelogged_yr0, twice_logged_delayed_fun)

twice_logged_delays <- do.call(rbind, predicted2L_values_list) %>% 
  cbind((harvest_delay = rep(harvest_delay_vt, each =rep_delay))) %>% rename(harvest_delay =3) %>% 
  cbind(original_habitat = "once-logged", 
        habitat = "twice-logged")  %>% 
  left_join(error_2L, by = "true_year") %>%  
  mutate(lwr_ACD = ACD - error95, 
         upr_ACD = ACD + error95) %>% 
  
  # #rough method of calculating ACD after second harvest assumes 0 ACD upon immediate reharest (e.g. with no time delay since first logging)
  # # this is unrealistic. Manually make it so that second harvest never leads to a lower ACD or error then 10 
  # # # this is wrong. But we also know from Riutta et al. that twice-logged forests are net sinks for first 10 years, so it's not mad to limit early year ACD recovery to artifially reduce carbon benefits
  #  mutate(
  #   lwr_ACD = ifelse(lwr_ACD < 10, 10, lwr_ACD),
  #   ACD = ifelse(ACD < 10, 10, ACD),
  #   upr_ACD = ifelse(upr_ACD < 10, 10, upr_ACD)) %>%
  
  #MANUALLY CAUSE CARBON TO PLATEU IN  FOREST ONCE IT REACHES OG VALUES 
  mutate(
    ACD = ifelse(ACD >max_val , NA, ACD),
    lwr_ACD = ifelse(ACD > max_val, NA, lwr_ACD),
    upr_ACD = ifelse(ACD > max_val, NA, upr_ACD)) %>%  
  #take on the max value if NA
  mutate(
    ACD = ifelse(is.na(ACD), max(ACD, na.rm = TRUE), ACD),
    lwr_ACD = ifelse(is.na(lwr_ACD), max(lwr_ACD, na.rm = TRUE), lwr_ACD),
    upr_ACD = ifelse(is.na(upr_ACD), max(upr_ACD, na.rm = TRUE), upr_ACD),
    time_since_intervention = true_year) %>% 
  
  #make true year correct 
  mutate(true_year = true_year + harvest_delay, 
         functional_habitat = "twice-logged") %>%  
  select(-error95)



#add time delays #### 

#NB - for forest that starts off as once-logged (i.e. where first harvesting happens in yr t-1, we are 
#only allowed to second harvest in yr 15). Yr 0 of 1L->2L in our scenarios is therefore actually 
#,technically, 15 yr once-logged forest.

#get once-logged -> once logged for the forest during the delay, before harvest 
delay_years_1L <- list()

for (i in harvest_delay_vt) {
  # Filter the dataframe based on the "true_age" value between 0 and i
  filtered_df <- oncelogged %>%
    filter(true_age >= 0, true_age <= i) %>%
    mutate(harvest_delay = i)
  
  # Add the filtered dataframe to the named list with a descriptive name
  delay_years_1L[[paste0( i)]] <- filtered_df
}

delay_years_1L <- bind_rows(delay_years_1L) %>% 
  # Remove rows where time_since_intervention is equal to harvest_delay
  filter(time_since_intervention != harvest_delay) %>%  
  mutate(habitat = "twice-logged") %>% 
  cbind(original_habitat = "once-logged") %>% 
  cbind(functional_habitat = "once-logged") %>%  
  rename(true_year= true_age)


#-----output of 3 ------
oncelogged_to_twice_logged <- rbind(twice_logged_delays,delay_years_1L) %>% 
  rename(functionalhabAge = time_since_intervention) %>% 
  filter(true_year < 61) %>%  
  group_by(harvest_delay) %>% 
  arrange(true_year) %>% ungroup

rownames(oncelogged_to_twice_logged) <- NULL
#---- # View each twice-logged 60yr typology schedule -------------------------


#---- Plotting 1&2 -----
p <- ggplot(from2L_2L_AND_Primary_2L, aes(x = true_year, y = ACD, color = original_habitat, linetype = original_habitat)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 250)) +
  theme_bw()

p

# Filter the data to get the label positions for "primary" habitat with "logged once here" and "logged twice here"
label_data_once <- from2L_2L_AND_Primary_2L %>%
  filter(original_habitat == "primary") %>%
  slice_min(true_year)

label_data_twice <- from2L_2L_AND_Primary_2L %>%
  filter(original_habitat == "twice-logged") %>%
  slice_min(true_year)

p + 
  geom_text(data = label_data_once, aes(label = "ReLogged after 30 yrs, original habitat is primary"), hjust = 0, vjust = 0) +
  geom_text(data = label_data_twice, aes(label = "Logged twice before start of scenario, oringal habitat is twice-logged"), hjust = 0, vjust = 1, color = "red")

#Looks good; minimum 15 years delay, then harvest, then
oncelogged_to_twice_logged  %>% 
  group_by(harvest_delay) %>% 
  ggplot(aes(x = true_year, y = ACD, colour = as.factor(harvest_delay))) +
  geom_line() +
  scale_y_continuous(limits = c(0, 250)) +
  theme_bw()

#------ print all twice-logged typologies ---------
oncelogged_to_twice_logged
from2L_2L_AND_Primary_2L


#--- Add deforested ----- 
deforested_c <- data.frame(true_age = seq(0,60), 
                           ACD = 0, 
                           lwr_ACD =0, 
                           upr_ACD = 0, 
                           habitat = "deforested") %>%  
  mutate(time_since_intervention = true_age)


#------ make all hab transitions ----------------
names(plantation)
names(primary_Vals)
names(L1_R)


combined <- plantation %>%  rbind(primary_Vals) %>% rbind(L1_R) %>% rbind(deforested_c) %>% rename(functional_habitat = habitat,
                                                                                                   functionalhabAge = time_since_intervention,
                                                                                                   true_year = true_age) 

#filter hab_years to remove the special cases we have manually calculated above 
hab_by_yearFilt <- hab_by_year  %>% filter(
  !(original_habitat == "once-logged" & habitat == "twice-logged") &
    !(original_habitat == "primary" & habitat == "twice-logged") &
    !(original_habitat == "twice-logged" & habitat == "twice-logged") #&
  # !(original_habitat == "once-logged" & habitat == "once-logged")
) 

hab_by_yearFilt %>% select(original_habitat,habitat) %>% unique

#stays twice logged
stays2L <- from2L_2L_AND_Primary_2L %>% filter(habitat == original_habitat) %>%  
  select(true_year, functional_habitat, ACD, lwr_ACD, upr_ACD) %>% mutate(functionalhabAge = true_year)

combined <- combined %>% rbind(stays2L) %>% 
  left_join(hab_by_yearFilt, by = c(  "functional_habitat" , "functionalhabAge", "true_year"), 
            relationship = "many-to-many")


#check correct number of rows for each hab transitions 
#XX <- combined %>% group_by(original_habitat, habitat) %>% count

# already has hab transitions 
oncelogged_to_twice_logged #already has harvest delays (only allows delays of >14...have to wait 15 years from first harvest in yr t-1)
from2L_2L_AND_Primary_2L #no harvest delays 


# ------ add harvest time-delay to scenarios -----------------
# define the delay before which habitat transition happens 
# here we apply a time window of 30 years, apply harvesting annually
delay <- seq(1:29)

#----create harvest delay function ----
#hab_by year currently assumes all harvesting starts in year 0. 
#what if delay harvests/first plantation plants? 
output_list <- vector("list", length(delay))

add_delay <- function(X) {
  
  # habTrans <-combined 
  #habTrans <- from2L_2L_AND_Primary_2L
  
  habTrans <- X  
  
  
  #filter data where habitat stays the same (e.g. once-logged stays once-logged)
  delay_with_zero <- c(0, delay)
  
  #store values for delay years, to join in for showing recovery in delay period 
  untransitioned_habs <- habTrans %>% filter(habitat == original_habitat) %>%  
    select(true_year, functional_habitat, ACD, lwr_ACD, upr_ACD) 
  
  #store habitats that undergo no transition 
  untransitioned_habs_delay <- habTrans %>% filter(habitat == original_habitat) %>%  
    crossing(delay_with_zero) %>% 
    rename(harvest_delay = delay_with_zero)
  
  #filter only data where there are transiions
  habTrans <- habTrans %>% filter(!habitat == original_habitat)
  
  hab_by_year0 <- habTrans %>% cbind(harvest_delay = 0)
  #ACD in original 
  for (i in delay) {
    #make a datafrmae of delayed years to add to beginning, to cause delay and leave habitat as original cover for the time of delay
    yearZ <-   habTrans %>%
      filter(true_year == 0) %>%
      uncount(i) %>%
      group_by(original_habitat,habitat) %>%
      mutate(true_year = seq_along(true_year)-1,
             functionalhabAge = seq_along(true_year)-1) %>%
      ungroup() %>% 
      #givve ACD in delay years an NA values   
      mutate(ACD = NA, 
             lwr_ACD = NA,
             upr_ACD = NA) %>%  
      
      #if ACD is NA, give the correct recovery delay years
      left_join(untransitioned_habs, by = c("true_year", "functional_habitat")) %>%
      mutate(
        ACD = ifelse(is.na(ACD.x), ACD.y, ACD.x),
        lwr_ACD = ifelse(is.na(lwr_ACD.x), lwr_ACD.y, lwr_ACD.x),
        upr_ACD = ifelse(is.na(upr_ACD.x), upr_ACD.y, upr_ACD.x)
      ) %>%
      select(-ACD.x, -lwr_ACD.x, -upr_ACD.x, -ACD.y, -lwr_ACD.y, -upr_ACD.y) %>% 
      filter(true_year <61) %>%  
      
      #catch cases where original_habitat is primary, and delay year carbon remains NA (caused by untransitioned_habs not having primary - primary transitioned)
      
      mutate(
        ACD = ifelse(original_habitat == "primary" &
                       functional_habitat == "primary" &
                       is.na(ACD),
                     203, ACD),
        lwr_ACD = ifelse(original_habitat == "primary" &
                           functional_habitat == "primary" &
                           is.na(lwr_ACD),
                         160, lwr_ACD),
        upr_ACD = ifelse(original_habitat == "primary" &
                           functional_habitat == "primary" &
                           is.na(upr_ACD),
                         250, upr_ACD)
      )
    
    #give age of original habitat during harvest delay, except for primary and deforested
    # mutate(functionalhabAge = case_when(
    #   functional_habitat != "primary" & functional_habitat != "deforested" ~ true_year,
    #   TRUE ~ functionalhabAge)) 
    
    #push true years by the length of the delay (e.g. add in 0s) and remove true year >60 
    delayed_df <-habTrans %>%   
      mutate(true_year = true_year + i ) %>% 
      filter(true_year <61)
    
    #combine then remove beyond 60th years 
    output <- yearZ %>%  rbind(delayed_df) %>% cbind(harvest_delay = paste(i))
    
    output_list[[i]] <- output
  }
  # 
  # #make sure we have all the years with same amount of rows
  #x <-  rbindlist(output_list)
  # test <- x %>% group_by(true_year) %>% count()
  # 
  # #this hab_by_year now includes the temporal dimension assuming different delays until first harvest
  habTrans<- rbindlist(output_list) %>% rbind(hab_by_year0)
  # 
  # #combine the transitioning and non-transitioning data 
  habTrans <- rbind(habTrans,untransitioned_habs_delay)
  
}

#----apply harvest delay function ----
from2L_2L_AND_Primary_2L <- add_delay(from2L_2L_AND_Primary_2L)
combined <- add_delay(combined)
XX <- combined %>% 
  group_by(original_habitat,habitat,harvest_delay)  %>% count

#----combine all habitat transitionssingle dataframe ----
allHabCarbon_60yrACD_withDelays <- combined %>% rbind(from2L_2L_AND_Primary_2L) %>% rbind(oncelogged_to_twice_logged)


#----final hard-coded correction----

# we are still missing all the delay year data from original twice-logged habitat
# mannually add this in
allHabCarbon_60yrACD_withDelays <- allHabCarbon_60yrACD_withDelays %>% left_join(stays2L, by = c("true_year", "functional_habitat","functionalhabAge")) %>%
  mutate(
    ACD = ifelse(original_habitat == "twice-logged" &
                   functional_habitat == "twice-logged" &
                   is.na(ACD.x),
                 ACD.y, ACD.x),
    lwr_ACD = ifelse(original_habitat == "twice-logged" &
                       functional_habitat == "twice-logged" &
                       is.na(lwr_ACD.x),
                     lwr_ACD.y, lwr_ACD.x),
    upr_ACD = ifelse(original_habitat == "twice-logged" &
                       functional_habitat == "twice-logged" &
                       is.na(upr_ACD.x),
                     upr_ACD.y, upr_ACD.x)) %>%
  select(-ACD.x, -lwr_ACD.x, -upr_ACD.x, -ACD.y, -lwr_ACD.y, -upr_ACD.y) %>%
  filter(true_year <61) %>% unique()

#check correct years for each habitat transition
XX <- allHabCarbon_60yrACD_withDelays %>% 
  group_by(original_habitat,habitat,harvest_delay)  %>% count


#-----add belowground carbon change----

#NEED TO CORRECT TO INCLUDE UNCERTAINTY CALCULATIONS; ATM ONLY CALCULATING THE MAIN BIT (NO UPR OR LWR)

#Assumptions made for incorporating belowground losses ####
#1. Adding belowground carbon and necromass. 
#2. For first 10 years, above ground ACD is offset by belowground losses 
#3. Belowground carbon recovers at the same rate as aboveground carbon thereafter (i.e. we  
#assume that at year 10, belowground suddenly becomes a sink, equivalent to aboveground) 
#4. Plantations lose all belowground carbon when deforested and don’t ever recover any belowground carbon. 
#5. Establishement of plantations on deforested ground doesn’t increase or decrease belowground carbon.  


#More formally
#conversion to plantation = JUST ACD (i.e. 0 belowground carbon)
#ACDforest + BCforest(t < 10) = ACDforest_t1 ----------> #  #if there is a habitat transition then first 10 years have the same ACD  ACD recovery is offset by belowground losses, or more formally
#ACDforest + BCforest(t > 10) = ACDforest_t + (ACDforest_t * 0.31)
#################
# #to test function below and see what's happening for replicatability 
# replica <- allHabCarbon_60yrACD_withDelays %>%  filter(harvest_delay == 15)
# x <- replica

belowground_fun <- function(x) {
  
  # 1. FOR PLANTATIONS 
  #for plantations, they lose all belowground carbon when deforested and don’t ever recover any belowground carbon. 
  #Establishement of plantations on deforested ground doesn’t increase or decrease belowground carbon.  
  #This is operationalised by belowground carbon always being 0 if functional habitat = plantation, so ACD is basically all that matters
  plantations <- x %>% 
    filter(str_detect(habitat, "albizia|eucalyptus")) %>% 
    mutate(full_carbon = ACD, 
           full_carbon_lwr = lwr_ACD, 
           full_carbon_upr = upr_ACD)
    
  #2. FOR ALL OTHER DATA 
 
   x <- x %>%  
    #remove plantation cases 
    filter(!str_detect(habitat, "albizia|eucalyptus"))

    #Get year 1 ACD 
  yr1_filt <- x %>% filter(functionalhabAge == 1) %>% 
    select(functional_habitat,original_habitat,habitat,harvest_delay,ACD, upr_ACD, lwr_ACD) %>%  
    rename(ACDt1 = ACD, 
           uprACDt1 = upr_ACD, 
           lwrACDt1 = lwr_ACD) %>% unique
  
  #if there is a habitat transition then first 10 years have the same fixed ACD 
  #as ACD recovery is offset by belowground losses, or more formally
  #ACDforest+BCforest(t < 10) = ACDforest_t1
  
  y <- x %>% left_join(yr1_filt) %>% 
      # Check conditions and mutate ACD accordingly
      mutate(
        full_carbon = case_when(
          
          #if there is a habitat transition and functional hab age <10, fixed full carbon as t1 ACD (equivalent to belowground carbon loss offsetting above ground gains)
          original_habitat != functional_habitat & functionalhabAge <= 10 ~ ACDt1,
          
          #if there is no habitat transition OR functional hab age >10, assume full carbon = ACD + BCD (where bcd = ACD*031)
          original_habitat == functional_habitat | functionalhabAge > 10 ~ ACD + (ACD* 0.31),
          TRUE ~ NA  # Otherwise NA
        ), 
        
        #get uppr and lwr bounds
        full_carbon_lwr = case_when(
          original_habitat != functional_habitat & functionalhabAge <= 10 ~ lwr_ACD,
                    original_habitat == functional_habitat | functionalhabAge > 10 ~ lwr_ACD + (lwr_ACD* 0.31),
          TRUE ~ NA), 
        
        full_carbon_upr = case_when(
          original_habitat != functional_habitat & functionalhabAge <= 10 ~ upr_ACD,
          original_habitat == functional_habitat | functionalhabAge > 10 ~ upr_ACD + (upr_ACD* 0.31),
          TRUE ~ NA)
      ) %>%  
    
    select(-c(ACDt1, uprACDt1, lwrACDt1))
  
  #recombine plantations and other data 
  full_data <- rbind(plantations, y)
  
}
  

allHabCarbon_60yrACD_withDelays <- belowground_fun(allHabCarbon_60yrACD_withDelays)

#----make master output ----

#reorder columns names 
names(allHabCarbon_60yrACD_withDelays)
allHabCarbon_60yrACD_withDelays <- allHabCarbon_60yrACD_withDelays %>%
  select(original_habitat, habitat, functional_habitat,
         true_year, functionalhabAge, harvest_delay,
         full_carbon, full_carbon_lwr, full_carbon_upr,
         ACD,lwr_ACD, upr_ACD)


#cursory plots ####
#nb only plot delay yrs 1 & 29 so that you can see what's going on
allHabCarbon_60yrACD_withDelays %>%
  filter(harvest_delay %in% c(1, 29)) %>% 
  ggplot(aes(x = true_year, y = full_carbon, colour = harvest_delay)) +
  geom_line() +
  facet_wrap(~original_habitat + habitat, scales = "free_y") +
  labs(x = "True Year", y = "Full Carbon") +
  theme_minimal()

#write Master output #####
#NB this outputs a master output where ACD refers only to the raw above ground carbon values and 
#where full_carbon incorporate belowground processes

write.csv(allHabCarbon_60yrACD_withDelays, "Outputs/allHabCarbon_60yr_withDelays.csv")

