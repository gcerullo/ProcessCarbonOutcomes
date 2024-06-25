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
library(cowplot)

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

#define params ####
harvest2ndACD_loss = 43.13944
#define ACD left after second rotation (Mg.ha-1), 30 years after 1st logging 
ACD_2L_yr0_30yrAfter1L_30yrAfter1L <-84.88416  # 30y3 1L =  128.02360   -  harvest2ndACD_loss
#define ACD after for scenarios beginning as 2L
ACD_2L_starting2L_15yrAfter1L <- 40.62932

#define inputs ####
#-----read in inputs showing habitat transitions --------
hab_by_year <- read.csv("Inputs/HabByYears.csv", strip.white = TRUE) %>%  
  rename(true_year = year, 
         functionalhabAge = functional_habAge, 
         habitat = transition_habitat) %>% select(-X) %>%  
  #remove all improved plantation yielding varieties (can add in if interested to)
  filter(!str_detect(habitat, "improved"))

  
#--------- read in plantation data ---------

plantation <- read.csv("Outputs/plantation_carbon.csv") %>% 
  select(-X) %>%  
  mutate(habitat = case_when(habitat == "albizia" ~ "albizia_current", 
                              habitat == "eucalyptus" ~"eucalyptus_current"))


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


# Plot 1L and Restored #### 
# add line indicating timing of restoration
# RestorLine <- data.frame(x=2, y=399, xend=23, yend=399, FACE="B: With active restoration")
# 
# plot1 <- ggplot(data= Logged, aes(YearsSinceLogging, ACD, FACE))+
#   xlab("Years Since Logging")+
#   ylab(expression(paste("Aboveground Carbon Density", " ","(Mg h", a^-1, sep = "", ")")))+
#   scale_y_continuous(limits=c(0,399))+
#   facet_grid(~FACE)+ # , scales="free_x", space="free"
#   # geom_point(col="grey25", pch="O")+ # grey30 grey20
#   geom_line(data= Logged,aes(YearsSinceLogging,ACD, group=Plot),  size=0.5, colour= "grey")+ # lines per plot
#   geom_line(data=NEWDAT,aes(YearsSinceLogging, ACD, group=FACE),  linetype=1, size=1.2, colour= "blue")+
#   geom_line(data=NEWDAT,aes(YearsSinceLogging,plo, group=FACE), size=1, linetype =3, colour= "blue")+
#   geom_line(data=NEWDAT,aes(YearsSinceLogging,phi, group=FACE),size=1, linetype =3,  colour= "blue")+
#   theme_bw()+ theme(strip.text.y = element_text( size = 12, hjust = .5), strip.text.x = element_text( size = 12, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
#   theme(strip.text.x = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
#   theme(strip.text.y = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
#   theme(axis.title.y = element_text(size = 14, angle = 90))+
#   theme(axis.title.x = element_text(size = 14, angle = 00))+
#   theme(axis.text = element_text(size = 14, colour = "black"))+labs(title="")+
#   geom_segment(data= RestorLine, mapping=aes(x=x, y=y, xend=xend, yend=yend),inherit.aes=FALSE, arrow=arrow(angle=65, ends="both", length=unit(0.1, "inches")), size=1, color="darkgreen")
# 
# 
# 
# #### Single panel version (easier to compare slopes)
# plot2 <- ggplot(data= Logged, aes(YearsSinceLogging, ACD, FACE))+
#   xlab("Years Since Logging")+
#   ylab(expression(paste("Aboveground Carbon Density", " ","(Mg h", a^-1, sep = "", ")")))+
#   scale_y_continuous(limits=c(0,399))+
#   geom_line(data= Logged,aes(YearsSinceLogging,ACD, group=Plot, colour =FACE),  size=0.5, alpha=0.2)+ # lines per plot
#   geom_line(data=NEWDAT,aes(YearsSinceLogging, ACD, colour =FACE),  linetype=1, size=1.5)+
#   geom_line(data=NEWDAT,aes(YearsSinceLogging,plo, colour =FACE), size=1, linetype =2)+
#   geom_line(data=NEWDAT,aes(YearsSinceLogging,phi, colour =FACE),size=1, linetype =2)+
#   theme_bw()+ theme(strip.text.y = element_text( size = 12, hjust = .5), strip.text.x = element_text( size = 12, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
#   theme(strip.text.x = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
#   theme(strip.text.y = element_text(size = 14, hjust = .5),strip.background=element_rect(colour="black", fill="white"))+
#   theme(axis.title.y = element_text(size = 14, angle = 90))+
#   theme(axis.title.x = element_text(size = 14, angle = 00))+
#   theme(axis.text = element_text(size = 14, colour = "black"))+labs(title="")+
#   scale_color_manual(values=c("blue", "red"))+ 
#   theme(legend.position = "none")



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

#get first 3yrs of 1L 
oncelogged <- L1_R %>% filter(habitat == "once-logged") %>% filter(true_age < 31)

#get restored and once logged properly organised 
L1_R <- L1_R %>% 
  select(-time_since_intervention) %>%  
  rename(functionalhabAge = true_age, 
         functional_habitat = habitat) 

#organise P - 1L and 1L - 1L transitions 
OnceLogged <- L1_R %>% filter(functional_habitat == "once-logged")
P_1L_df <- OnceLogged %>% mutate(original_habitat = "primary", 
                                 habitat = "once-logged")

L1_1L_df <- OnceLogged %>% mutate(original_habitat = "once-logged", 
                                  habitat = "once-logged")

#organise P - R and 1L-R transitions 
Restored <- L1_R %>% filter(functional_habitat == "restored") 
P_R_df <- Restored %>% mutate(original_habitat = "primary", 
                              habitat = "restored")
L1_R_df <- Restored %>% mutate(original_habitat = "once-logged", 
                              habitat = "restored")  

#------------- get correct years -------------
#make corrections so that we have values for all (and for the correct) years


#----60 yrs of plantations ----

#1.make plantation carbon over 60 years
plantation <-  plantation %>% 
  rename(functionalhabAge = plantationAge)  

#add in zero plantation age (for the very first year of clearance)
zero_eucPlant <- data.frame(ACD = 0, 
                            lwr_ACD = 0, 
                            upr_ACD = 0, 
                            habitat = "eucalyptus_current", 
                            functionalhabAge = 0)
zero_albPlant <- data.frame(ACD = 0, 
                            lwr_ACD = 0, 
                            upr_ACD = 0, 
                            habitat = "albizia_current", 
                            functionalhabAge = 0)


plantation_df <- zero_eucPlant %>% rbind(zero_albPlant) %>% 
  rbind(plantation) %>% 
  rename(functional_habitat = habitat)

#select all hab crosses possible for plantation 
habcrossAlb <- hab_by_year %>% 
  filter(functional_habitat == "albizia_current") %>%  
  select(original_habitat, habitat) %>% 
  unique()

#select all hab crosses possible for plantation 
habcrossPlantEuc <- hab_by_year %>% 
  filter(functional_habitat =="eucalyptus_current") %>%  
  select(original_habitat, habitat) %>% 
  unique()


#expand grid so that we have functionalhabAge for all transitions
alb_df <- plantation_df %>% 
  filter(functional_habitat =="albizia_current") %>%
  expand_grid(habcrossAlb)

euc_df <- plantation_df %>% 
  filter(functional_habitat =="eucalyptus_current") %>%
  expand_grid(habcrossPlantEuc)

plant_df <- rbind(alb_df,euc_df)


#----60 yrs of primary ----
#2.make primary carbon over 60 years

primary <- primary_Vals %>% 
  rename(functional_habitat = habitat) %>%  
  cbind(functionalhabAge = 0)

primary_cross <- hab_by_year %>% 
  filter(functional_habitat == "primary") %>%  
  select(original_habitat, habitat) %>% 
  unique()

primary_df <- expand_grid(primary,primary_cross)

#--- Add deforested ----- 
deforested <- data.frame(functionalhabAge = 0, 
                           ACD = 0, 
                           lwr_ACD =0, 
                           upr_ACD = 0, 
                           functional_habitat = "deforested")


deforested_cross <- hab_by_year %>% 
  filter(functional_habitat == "deforested") %>%  
  select(original_habitat, habitat) %>% 
  unique()

deforested_df <- expand_grid(deforested,deforested_cross)

#combine all carbon except for cases with more complex 2L typologies, defined below
carbhabs <- deforested_df %>% 
  rbind(primary_df) %>%
  rbind(plant_df) %>%
  rbind(P_R_df) %>% 
  rbind(L1_R_df) %>%
  rbind(P_1L_df) %>% 
  rbind(L1_1L_df)


#----60 yrs of twice-logged ----
#--- three 2L typologies (1,2,3) ----

#Briefly: 
#1 = P -> 2L
#2 = 2L -> 2L 
#3 = 1L -> 2L 

#NOTE THERE ARE THREE TYPES OF TWICE-LOGGING CURVES THRU TIME BUILT HERE: 
#1. ACD_2L_yr0_30yrAfter1L_30yrAfter1L - second rotations happens 30 years after 1st logging. Always the case in all_primary and all_primary_deforested scenarios

#2, ACD_2L_starting2L_15yrAfter1L - 
#here the second rotation happens 15 years after 1st rotation 
#( the second rotation is assumed to have happened in t-1, before our scenarios start) 
# We incoroprate this for all scenarios that start 2L - as this more rapid re-harvest of forest more closely resembles the true 
#logging dynamics in the landscape where we surveyed twice logged forests.


#i.e in   t-16----t-15 --------------------------t-1--------t-0------------------------------------t-60
#         P        1L                             2L        2L                                      2L

#We also make this true for scenarios that go from 1L (at the beginning of the scenarios) -> 2L....
#....which is operationalised by not allowing 1L-2L transitions for 15 years....see below

#3. twice_logged_delays  - second rotation happens at varying intervals after the once-logging from 15-30 years, depending on the harvesting delay. 
# this is important for scenarios that begin with once-logged forest (eg. mostly_1L and mostly_1L_deforested scenarios) 
# the delay matters because if the second harvest happens after a big time delay (e.g. 25 yrs) then ACD was recovering all this time in once-logged, so the 2nd harvest leaves behind higher ACD 
#NB - no logging of 1L -> 2L is allowed in the first 15 yrs, as the forest was logged in t-1. 

#define some key params and functions ####

#Above ground carbon in twice-logged straight after logging (y-intercept)
print(ACD_2L_yr0_30yrAfter1L_30yrAfter1L)  #this assumes second logging happens 30 years after once-logging 
print(ACD_2L_starting2L_15yrAfter1L) #this assume second logging happens 15 yrs after first logging, at year t-1...this is the ACD for scenarios starting as twice logged

# Define the y-intercept, slope, and confidence intervals
Slope_1L_df <- L1_R %>% filter(functional_habitat == "once-logged")   
Slope_1L <-  lm(ACD ~ functionalhabAge, data = Slope_1L_df)
# Extract the slope using tidy()
Slope_1L <- tidy(Slope_1L)$estimate[2]  
years <- 0:60

# Calculate the predicted y-values for each year
#1.
predicted2L_values <- ACD_2L_yr0_30yrAfter1L_30yrAfter1L + Slope_1L * years # #this assumes second logging happens 30 years after once-logging

#2.
predicted2L_values_starting <- ACD_2L_starting2L_15yrAfter1L + Slope_1L * years#this assume second logging happens 15 yrs after first logging, at year t-1

# Calculate the lower and upper limits of the confidence intervals
#assume sameCI as once-logged
error_2L <- L1_R %>% filter(functional_habitat == "once-logged") %>%
  mutate(error95 = ACD - lwr_ACD) %>%
  select(error95) %>%  
  cbind(true_year = 0:60)

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

## for 1: P -> 2L ####
twice_L_30yrPost <-  twice_log_fun(predicted2L_values) %>% 
  cbind(original_habitat = "primary")%>%
  select(-true_age)

#adjust primary -> 2L so that is undergoes a once-logged conversion, then a twice-logged conversion after 30 yrs 
names(twice_L_30yrPost)

# yr0_primary_2L <- primary_Vals  %>% slice(1) %>% 
# #rename(true_year = true_age)  %>% 
#   mutate(habitat = "twice-logged",
#          original_habitat = "primary", 
#          functional_habitat = "primary",
#          functionalhabAge = 0) 


yr1_30_primary_2L <- oncelogged  %>% rename(true_year = true_age) %>%
  filter(true_year >= 1, true_year <= 29) %>%  
  mutate(habitat = "twice-logged",
         original_habitat = "primary", 
         functional_habitat = "once-logged") %>%  
  rename(functionalhabAge = true_year) %>%  
  select(-time_since_intervention)

# beginning_habs_primary_2L <- rbind(yr0_primary_2L,yr1_30_primary_2L) 

twice_L_30yrPost <- twice_L_30yrPost %>% mutate(true_year = true_year + 30) %>%
  rename(functionalhabAge = time_since_intervention) %>% 
  filter(true_year <61) %>% 
  mutate(functional_habitat = "twice-logged") %>% 
  select(-true_year)

#final P_2L (i.e. 1. df)
P_2L_df <- rbind(yr1_30_primary_2L,twice_L_30yrPost) 

##for 2.(2L -> 2L) #### 
L2_2L_df <-  twice_log_fun(predicted2L_values_starting) %>% cbind(original_habitat = "twice-logged") %>%  
  cbind(functional_habitat = "twice-logged") %>% 
  rename(functionalhabAge = time_since_intervention) %>% 
  select(-c(true_age, true_year))

##for 3. (1L->2L)####
#nb we can't harvest once-logged forest until year 15 in our scenarios
#because once-logging happeneng -t15 years before scenario start 
#for parcels beginning as once-logged 

oncelogged <- L1_R %>%
  filter(functional_habitat == "once-logged") %>% 
  filter(functionalhabAge >14) %>% 
  filter(functionalhabAge <45) %>% 
  #once-logged forest starting in our scenarios is actually 15 yrs old already 
  #so this corrects for this. 
  mutate(functionalhabAge = functionalhabAge - 15) %>% 
  mutate(original_habitat = "once-logged", 
         habitat = "twice-logged") 

twicelogged <- P_2L_df %>% 
  filter(functional_habitat == "twice-logged") %>% 
  mutate(original_habitat = "once-logged") 


#nb - remember, when operationalising this transition in
#next script, remember to filter out any delay <15
#this will ensure that we enable once-logged forest to regrow for a min
#of 15 from scenario start (= to 30 yrs since harvest, at t-15)
L1_L2_df <- rbind(oncelogged, twicelogged)
  
#---- # View each twice-logged 60yr typology schedule -------------------------


#---- Plotting 1&2 -----
p3 <- ggplot(P_2L_df , aes(x = functionalhabAge, y = ACD, color = functional_habitat, linetype = original_habitat)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 250)) +
  theme_bw()

p2<- ggplot(L2_2L_df , aes(x = functionalhabAge, y = ACD, color = functional_habitat, linetype = original_habitat)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 250)) +
  theme_bw()

p3<- ggplot(L1_L2_df , aes(x = functionalhabAge, y = ACD, color = functional_habitat, linetype = original_habitat)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 250)) +
  theme_bw()

plot_grid(p1, p2,p3, col =2) # looks good (remember for 1L-2L, the once logged forest starts in functionhabAge = 0 at actually 15 yrs old)

#------ combine all hab transitions ----------------
#24.06.24!!!!####
#COME BACK TO HERE ####
carbhabs <- carbhabs %>%  
  rbind(P_2L_df) %>% 
  rbind(L2_2L_df) %>% 
  rbind(L1_L2_df)


#check correct number of rows for each hab transitions 
XX <- carbhabs %>% group_by(original_habitat, habitat) %>% count

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

