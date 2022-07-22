# Load required packages

library(dplyr)
library(reshape2)

library(devtools)
install_github("dsrobertson/onlineFDR")
library(onlineFDR)

library(ggplot2)
library(ggnewscale)
library(patchwork)

################################################################################
############################   Test levels    ##################################
################################################################################

set.seed(3)

N = 300

F1 = rnorm(N, mean = 3, sd = 1)

pi.theta = rbinom(N, 1, 0.5)

theta = rep(0,N)
theta[pi.theta==1] = F1[pi.theta==1]

Z = theta + rnorm(N)

pval = pnorm(-Z)


# LORD
alphai.LORD = log10(LORD(pval)$alphai)

# SAFFRON
alphai.SAFFRON = log10(SAFFRON(pval)$alphai)

# ADDIS
alphai.ADDIS = log10(ADDIS(pval)$alphai)


# Alpha-spending
alphai.Bonf = log10(Alpha_spending(pval)$alphai)


# Create plot

alphai = data.frame(Hypothesis = 1:N,
                    LORD = alphai.LORD, SAFFRON = alphai.SAFFRON,
                    ADDIS = alphai.ADDIS)

alphai = melt(alphai, id.vars = 'Hypothesis')


cbbPalette = c("#E69F00", "#56B4E9", "#009E73")

alphai.plot = ggplot(data = alphai, aes(x = Hypothesis, y = value)) +
  geom_line(aes(colour = variable), show.legend = FALSE) + 
  facet_wrap(.~variable) + 
  scale_colour_manual(values=cbbPalette) + 
  new_scale_color() + 
  geom_hline(data = data.frame(yint = log10(0.05),
                               variable = factor(c('LORD', 'SAFFRON', 'ADDIS'),
                                                 levels = levels(alphai$variable))),
             aes(yintercept = yint, linetype = 'dashed'), size = 1) + 
  geom_line(data = data.frame(Hypothesis = 1:N,
                              value = alphai.Bonf,
                              variable = factor(c('LORD', 'SAFFRON', 'ADDIS'),
                                                levels = levels(alphai$variable))),
            aes(x = Hypothesis, y = value, colour = 'Bonf'),
            linetype = 'dashed', size = 1) +
  scale_colour_manual(values="#CC79A7",
                      labels = 'Alpha-spending') +
  xlab('t') + 
  scale_x_continuous(limits = c(1,300)) +
  theme(legend.key.size=unit(1.3,"cm")) +
  ylab(expression(paste(log[10],' (',alpha[t], ')'))) + labs(color='') + 
  scale_linetype_manual(name = "", values = 2, labels = 'Uncorrected') +
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 0)) + 
  theme(text = element_text(size = 16), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


ggsave("alphai_plot.png", alphai.plot, width = 20, height = 7, units = "cm")


################################################################################
############################ Simulation study ##################################
################################################################################

set.seed(3)

nrep = 10^4

pi.vec = c(seq(0.01, 0.09, by = 0.01), seq(0.1, 0.9, by = 0.05))

bonf.power = bonf.FDR = BH.power = BH.FDR = LORD.power = LORD.FDR = 
  SAFFRON.power = SAFFRON.FDR = ADDIS.power = ADDIS.FDR = 
  Uncorrected.power = Uncorrected.FDR = rep(0, length(pi.vec))

bonf.TP = bonf.FDP = BH.TP = BH.FDP = LORD.TP = LORD.FDP = 
  SAFFRON.TP = SAFFRON.FDP = ADDIS.TP = ADDIS.FDP = 
  Uncorrected.TP = Uncorrected.FDP = rep(0, nrep)

################################################################################

N = 1000


for (i in 1:length(pi.vec)){
  
  for(j in 1:nrep) {
    
    F0 = rnorm(N, mean = -0.5, sd = 0.1)
    F1 = rnorm(N, mean = 3, sd = 1)
    
    pi.theta = rbinom(N, 1, pi.vec[i])
    
    theta = F0
    theta[pi.theta==1] = F1[pi.theta==1]
    
    Z = theta + rnorm(N)
    
    pval = pnorm(-Z)
    
    
    # Alpha-spending
    R.bonf = Alpha_spending(pval)$R
    
    bonf.TP[j] = sum(R.bonf[theta > 0] == 1)/sum(theta > 0)
    
    if(sum(R.bonf == 1) > 0) {
      bonf.FDP[j] = sum(R.bonf[theta < 0] == 1)/sum(R.bonf == 1)
    } else { bonf.FDP[j] = 0 }
    
    
    # BH
    R.BH = (p.adjust(pval, method = "BH") <= 0.05)
    
    BH.TP[j] = sum(R.BH[theta > 0] == 1)/sum(theta > 0)
    
    if(sum(R.BH == 1) > 0) {
      BH.FDP[j] = sum(R.BH[theta < 0] == 1)/sum(R.BH == 1)
    } else { BH.FDP[j] = 0 }
    
    
    # LORD
    R.LORD = LORD(pval)$R
    
    LORD.TP[j] = sum(R.LORD[theta > 0] == 1)/sum(theta > 0)
    
    if(sum(R.LORD == 1) > 0) {
      LORD.FDP[j] = sum(R.LORD[theta < 0] == 1)/sum(R.LORD == 1)
    } else { LORD.FDP[j] = 0 }
    
    
    # SAFFRON
    R.SAFFRON = SAFFRON(pval)$R
    
    SAFFRON.TP[j] = sum(R.SAFFRON[theta > 0] == 1)/sum(theta > 0)
    
    if(sum(R.SAFFRON == 1) > 0) {
      SAFFRON.FDP[j] = sum(R.SAFFRON[theta < 0] == 1)/sum(R.SAFFRON == 1)
    } else { SAFFRON.FDP[j] = 0 }
    
    
    # ADDIS
    R.ADDIS = ADDIS(pval)$R
    
    ADDIS.TP[j] = sum(R.ADDIS[theta > 0] == 1)/sum(theta > 0)
    
    if(sum(R.ADDIS == 1) > 0) {
      ADDIS.FDP[j] = sum(R.ADDIS[theta < 0] == 1)/sum(R.ADDIS == 1)
    } else { ADDIS.FDP[j] = 0 }
    
    
    # Uncorrected
    R.Uncorrected = (pval <= 0.05)
    
    Uncorrected.TP[j] = sum(R.Uncorrected[theta >0] == 1)/sum(theta >0)
    
    if(sum(R.Uncorrected == 1) > 0) {
      Uncorrected.FDP[j] = sum(R.Uncorrected[theta < 0] == 1)/sum(R.Uncorrected == 1)
    } else { Uncorrected.FDP[j] = 0 }
    
  }
  
  bonf.power[i] = mean(bonf.TP)
  bonf.FDR[i] = mean(bonf.FDP)
  BH.power[i] = mean(BH.TP)
  BH.FDR[i] = mean(BH.FDP)
  LORD.power[i] = mean(LORD.TP)
  LORD.FDR[i] = mean(LORD.FDP)
  SAFFRON.power[i] = mean(SAFFRON.TP)
  SAFFRON.FDR[i] = mean(SAFFRON.FDP)
  ADDIS.power[i] = mean(ADDIS.TP)
  ADDIS.FDR[i] = mean(ADDIS.FDP)
  Uncorrected.power[i] = mean(Uncorrected.TP)
  Uncorrected.FDR[i] = mean(Uncorrected.FDP)
  
}


FDR = data.frame(pi.vec, Uncorrected.FDR, BH.FDR,
                 LORD.FDR, SAFFRON.FDR, ADDIS.FDR,
                 bonf.FDR)

FDR = melt(FDR, id.vars = 'pi.vec')

pow = data.frame(pi.vec, Uncorrected.power, BH.power, 
                 LORD.power, SAFFRON.power, ADDIS.power,
                 bonf.power)


pow = melt(pow, id.vars = 'pi.vec')

cbbPalette = c("#000000", "#0072B2", "#E69F00", "#56B4E9", "#009E73",
                "#CC79A7") #"#D55E00", "#A020F0"

leg_label = c("Uncorrected","BH", "LORD", "SAFFRON", "ADDIS",
              "Alpha-spending")

FDR.plot1000 = ggplot(data=FDR,
                     aes(x=pi.vec, y=value,
                         colour = variable, linetype = variable)) +
  geom_line(size = 1) + 
  scale_colour_manual(values=cbbPalette,
                      name = "          Algorithms",
                      labels = leg_label) +
  scale_linetype_manual(values=c(2, 2, 1, 1, 1, 2),
                        name = "          Algorithms",
                        labels = leg_label) +
  theme(legend.key.size=unit(1.3,"cm")) +
  xlab(expression(pi[1])) + ylab("FDR") +
  scale_x_continuous(limits = c(0,0.91),
                     expand = c(0, 0),
                     breaks = seq(0, 0.9, by = 0.1)) +
  scale_y_continuous(limits = c(0, 0.7),
                     expand = c(0, 0)) +
  geom_hline(yintercept=0.05, colour = 'red')

power.plot1000 = ggplot(data=pow,
                       aes(x=pi.vec, y=value,
                           colour = variable, linetype = variable)) +
  geom_line(size = 1) + 
  scale_colour_manual(values=cbbPalette,
                      name = "          Algorithms",
                      labels = leg_label) +
  scale_linetype_manual(values=c(2, 2, 1, 1, 1, 2),
                        name = "          Algorithms",
                        labels = leg_label) +
  theme(legend.key.size=unit(1.3,"cm")) +
  xlab(expression(pi[1])) + ylab("Power") +
  scale_x_continuous(limits = c(0,0.91), expand = c(0, 0),
                     breaks = seq(0, 0.9, by = 0.1))+
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     breaks = seq(0, 1, by = 0.1))


# Arrange plots

power.FDR.plot = (FDR.plot1000 / power.plot1000) + 
  plot_layout(nrow = 2, guides = "collect")


# Save plots

ggsave("power_plot.png", power.plot1000, width = 18, height = 12, units = "cm")
ggsave("FDR_plot.png", FDR.plot1000, width = 18, height = 12, units = "cm")
ggsave("power_FDR_plot.png", power.FDR.plot, width = 17, height = 20, units = "cm")
