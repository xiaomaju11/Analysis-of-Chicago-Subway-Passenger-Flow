#############################################
# 0. Nettoyage de l'environnement
#############################################
rm(list = ls())

# Chargement des bibliothèques nécessaires
library(zoo)
library(timeDate)
library(forecast)
library(xts)
library(mgcv)
library(tidyverse)
library(splines)

#############################################
# 1. Lecture et préparation des données
#############################################
# Lecture du fichier CSV dont les colonnes sont :
# 1) service_date
# 2) day_type
# 3) bus
# 4) rail_boardings
# 5) total_rides
data <- read.csv("D:/project-R/Données_qualité_air_fréquentation/CTA_-_Ridership_-_Daily_Boarding_Totals.csv", 
                 header = TRUE, stringsAsFactors = FALSE)

# Vérifier les noms de colonnes
print(names(data))
# On suppose qu'elles correspondent exactement à :
# [1] "service_date"   "day_type"   "bus"   "rail_boardings"   "total_rides"

# Conversion de la date (format supposé "dd/mm/yyyy")
data$service_date <- as.Date(data$service_date, format = "%d/%m/%Y")

# Filtrer uniquement les données avant 2020-01-01
data <- subset(data, service_date < as.Date("2020-01-01"))
summary(data)

#############################################
# 2. Agrégation mensuelle (uniquement avant 2020)
#############################################
# Créer la variable "YearMonth" pour agréger par mois
data$YearMonth <- format(data$service_date, "%Y-%m")

# Calcul de la moyenne mensuelle de total_rides
monthly_avg <- aggregate(total_rides ~ YearMonth, data = data, FUN = mean)

# Conversion de YearMonth en date (on prend le 1er jour du mois)
monthly_avg$Date <- as.Date(paste0(monthly_avg$YearMonth, "-01"))
# Tri par date
monthly_avg <- monthly_avg[order(monthly_avg$Date), ]

# Construction d'un objet ts (fréquence = 12) et d'un objet xts
start_year <- as.numeric(format(min(monthly_avg$Date), "%Y"))
start_mon  <- as.numeric(format(min(monthly_avg$Date), "%m"))
monthly_ts <- ts(monthly_avg$total_rides, start = c(start_year, start_mon), frequency = 12)
monthly_xts <- xts(monthly_avg$total_rides, order.by = monthly_avg$Date)
plot(monthly_xts, main= "Valeur moyenne")


# Visualisation : série mensuelle, boxplot et histogramme
par(mfrow = c(1, 2))
plot(monthly_ts, xlab = "Temps", ylab = "Total Rides (avant 2020)",
     main = "Série mensuelle (avant 2020)")
boxplot(as.numeric(monthly_ts), main = "Boxplot mensuel", ylab = "Total Rides")
hist(monthly_ts, main = "Histogramme mensuel", xlab = "Total Rides")
par(mfrow = c(1, 1))

#############################################
# 3. Analyse de la tendance (fenêtre = 12)
#############################################

# On prend l'ensemble des données mensuelles
n <- length(monthly_ts)
t <- 1:n

## 3.1 Régression non linéaire (différents exposants de t)
reg <- lm(monthly_ts ~ t+I(t^2) + I(t^3))
trend_fit <- reg$fitted.values



# Visualisation des résidus et ACF
par(mfrow = c(1, 2))
plot(as.numeric(reg$residuals), type = "l", main = "Résidus de la régression", ylab="Résidus")
acf(as.numeric(reg$residuals), lag.max = 50, main = "ACF des résidus")

# Conversion en xts pour tracer sur le même axe de temps
trend_xts <- xts(as.numeric(trend_fit), order.by = monthly_avg$Date)
plot(monthly_xts, type = "l", main = "Tendance (régression non linéaire)",
     xlab = "Temps", ylab = "Total Rides")
lines(trend_xts, col = "red", lwd = 2)
print(summary(reg))

## 3.2 Tendance par moyenne mobile (fenêtre = 12)

ma_window <- 12
ma_trend <- stats::filter(monthly_ts, filter = rep(1/ma_window, ma_window),
                          method = "convolution", sides = 2, circular = FALSE)
ma_trend_xts <- xts(as.numeric(ma_trend), order.by = monthly_avg$Date)
plot(monthly_xts, type = "l", main = "Tendance (moyenne mobile, fenêtre=12)",
     xlab = "Temps", ylab = "Total Rides")
lines(ma_trend_xts, col = "blue", lwd = 2)

## 3.3 Tendance par noyau gaussien (h = 12)
h <- 12
x_seq <- seq(1, n, length.out = n)
gauss_kernel <- function(x_val) {
  dnorm(x_val - t, 0, sd = sqrt(h/2)) / sum(dnorm(x_val - t, 0, sd = sqrt(h/2)))
}
W <- matrix(unlist(lapply(x_seq, gauss_kernel)), nrow = n, ncol = n, byrow = FALSE)
kernel_trend <- colSums(as.numeric(monthly_ts) * W)
kernel_trend_xts <- xts(kernel_trend, order.by = monthly_avg$Date)
plot(monthly_xts, type = "l", main = "Tendance (noyau gaussien, h=12)",
     xlab = "Temps", ylab = "Total Rides")
lines(kernel_trend_xts, col = "red", lwd = 2)

## 3.4 Tendance par GAM (splines)
t_numeric <- as.numeric(monthly_avg$Date)
g <- gam(monthly_ts ~ t_numeric + s(t_numeric, k = 50))
gam_trend <- g$fitted.values
gam_trend_xts <- xts(as.numeric(gam_trend), order.by = monthly_avg$Date)
plot(monthly_xts, type = "l", main = "Tendance (GAM)", 
     xlab = "Temps", ylab = "Total Rides")
lines(gam_trend_xts, col = "red", lwd = 2)
print(summary(g))

#############################################
# 4. Analyse de la saisonnalité (fenêtre=12) & prévision
#############################################
# 4.1 Série détendrée (on enlève la tendance non linéaire par exemple)
detrended <- monthly_xts - gam_trend_xts
plot(detrended, type = "l", main = "Série détendrée (données - tendance)",
     xlab = "Temps", ylab = "Résidus")
acf(as.numeric(detrended), lag.max = 50, main = "ACF de la série détendrée")

# 4.2 Estimation de la saisonnalité (méthode 1 : moyenne mobile, fenêtre=12)
K <- 12
mb.season <- stats::filter(detrended, filter = rep(1/K, K), 
                           method = "convolution", sides = 2, circular = TRUE)
mb.season_xts <- xts(as.numeric(mb.season), order.by = monthly_avg$Date)
plot(detrended, type = "l", main = "Saisonnalité (moyenne mobile, K=12)")
lines(mb.season_xts, col = "blue", lwd = 2)

# 4.3 Estimation de la saisonnalité (méthode 2 : noyau tricube, fenêtre=12)
h_season <- 12
x_seq2 <- seq(1, n, length.out = n)
tricube_kernel <- function(u) {
  w <- (1 - abs(u)^3)^3
  w[abs(u) >= 1] <- 0
  return(w)
}
W_season <- matrix(unlist(lapply(x_seq2, function(x_val) {
  u <- (x_val - t) / sqrt(h_season/2)
  w <- tricube_kernel(u)
  w / sum(w)
})), ncol = n, nrow = n, byrow = FALSE)
season_kernel <- colSums(as.numeric(detrended) * W_season)
season_kernel_xts <- xts(season_kernel, order.by = monthly_avg$Date)
plot(detrended, type = "l", main = "Saisonnalité (noyau tricube, h=12)",
     xlab = "Temps", ylab = "Effet saisonnier")
lines(season_kernel_xts, col = "red", lwd = 2)
#4.4 Estimation de la saisonnalité (méthode 3 : Fourier, fenêtre=12)
w=2*pi/12##windows
fourier<-cbind(cos(w*t), sin(w*t))
K<-20
for(i in c(2:K))
{
  fourier<-cbind(fourier,cos(i*w*t), sin(i*w*t))
}
matplot(fourier[,c(1,3)],type='l')
dim(fourier)


r2<- NULL
err <- NULL
for(i in seq(2, ncol(fourier), by=2))
{
  reg<-lm(detrended~fourier[,1:i]-1)
  s <- summary(reg)
  r2  <- c(r2, s$r.squared)
  err <- c(err, sum(reg$residuals^2))
}

plot(err, type='b', pch=20, ylim=c(50, 200))
lines(err+2*c(1:ncol(fourier)), col='red')


reg<-lm(detrended~fourier[,1:20]-1)
summary(reg)

ychap.lm.season<-xts(as.numeric(reg$fitted),order.by=monthly_avg$Date)
plot(detrended,type='l')
lines(ychap.lm.season,col='red', lwd=2)
# Vérification de la série résiduelle
par(mfrow = c(1, 1))
resid_final <- detrended - ychap.lm.season
resid_final_filled <- na.approx(resid_final)
plot(resid_final_filled, type = "l", main = "Série résiduelle finale",
     xlab = "Temps", ylab = "Résidus")
acf(as.numeric(resid_final_filled), lag.max = 50, main = "ACF de la série résiduelle finale")

# (Optionnel) Deuxième estimation de la saisonnalité
# ... selon le même principe, si besoin

#############################################
# 5. STL 分解 & ARIMA 建模
#############################################
plot(monthly_ts, main = "Série mensuelle (avant 2020)", xlab = "Temps", ylab = "Total Rides")
stl_decomp <- stl(monthly_ts, s.window = "periodic")
plot(stl_decomp, main = "Décomposition STL")
trend_stl <- stl_decomp$time.series[, "trend"]
seasonal_stl <- stl_decomp$time.series[, "seasonal"]
remainder_stl <- stl_decomp$time.series[, "remainder"]
acf(remainder_stl, lag.max = 50, main = "ACF des résidus (STL)")
fit_arima <- auto.arima(remainder_stl)
summary(fit_arima)
checkresiduals(fit_arima)

#############################################
# 6. Lissage exponentiel et prévision
#############################################
# Fonctions de lissage exponentiel (restent inchangées)
expSmooth <- function(x, alpha) {
  xs <- numeric(length(x))
  xs[1] <- x[1]
  for(i in 2:length(x)) {
    xs[i] <- (1 - alpha) * xs[i-1] + alpha * x[i]
  }
  xs
}

DoubleExpSmooth <- function(x, alpha) {
  n <- length(x)
  l <- numeric(n)
  b <- numeric(n)
  smooth <- numeric(n)
  l[1] <- x[1]
  b[1] <- if(n > 1) x[2] - x[1] else 0
  smooth[1] <- x[1]
  for(i in 2:n) {
    l[i] <- smooth[i-1] + (1 - (1 - alpha)^2) * (x[i] - smooth[i-1])
    b[i] <- b[i-1] + alpha^2 * (x[i] - smooth[i-1])
    smooth[i] <- l[i] + b[i]
  }
  list(smooth = smooth, l = l, b = b)
}

alpha <- 0.5
# Calcul du lissage exponentiel simple
sm_simple <- expSmooth(monthly_avg$total_rides, alpha)
plot(monthly_avg$Date, monthly_avg$total_rides, type = "l", col = "red",
     xlab = "Temps", ylab = "Total Rides", main = "Lissage exponentiel simple")
lines(monthly_avg$Date, sm_simple, col = "blue")
legend("topleft", legend = c("Données", "Lissage simple"), col = c("red", "blue"), lty = 1)

# Calcul du lissage exponentiel double
sm_double <- DoubleExpSmooth(monthly_avg$total_rides, alpha)
plot(monthly_avg$Date, monthly_avg$total_rides, type = "l", col = "red",
     xlab = "Temps", ylab = "Total Rides", main = "Lissage exponentiel double")
lines(monthly_avg$Date, sm_double$smooth, col = "blue")
legend("topleft", legend = c("Données", "Lissage double"), col = c("red", "blue"), lty = 1)

# Affichage des paramètres (l et b)
plot(monthly_avg$Date, sm_double$l, type = "l", col = "red",
     ylim = range(c(sm_double$l, sm_double$b)),
     main = "Paramètres du lissage double")
lines(monthly_avg$Date, sm_double$b, col = "black")
legend("topleft", legend = c("l", "b"), col = c("red", "black"), lty = 1)

## 6.1 Prévision améliorée avec axe temporel
# Fonction de prévision (identique)
predictExpSmooth <- function(smooth, l, b, inst, horizon) {
  c(smooth[1:inst], l[inst] + b[inst] * (1:horizon))
}
# Définition de la taille de la période d'entraînement et de l'horizon de prévision
inst <- 228      # Utiliser les 80 premiers mois pour l'entraînement
horizon <- 12   # Prévoir les 12 mois suivants

# Calcul de la prévision
data_predict <- predictExpSmooth(sm_double$smooth, sm_double$l, sm_double$b,
                                 inst = inst, horizon = horizon)

# Construction des dates pour la prévision :
# - La partie historique : monthly_avg$Date[1:inst]
# - La partie prévision : construire une séquence à partir du mois (inst+1)
dates_hist <- monthly_avg$Date[1:inst]
dates_forecast <- seq(from = monthly_avg$Date[inst] + 1, by = "month", length.out = horizon)
# Tracer la prévision sur le même graphique en utilisant l'axe des dates
plot(monthly_avg$Date, monthly_avg$total_rides, pch = 20,
     xlab = "Date", ylab = "Total Rides", main = "Prévision (lissage double)",
     ylim = range(c(monthly_avg$total_rides, data_predict), na.rm = TRUE))
# Tracer la partie historique prédit (en rouge)
lines(monthly_avg$Date[1:inst], data_predict[1:inst], col = "red", lwd = 2)
# Tracer la partie prévision (en bleu)
lines(dates_forecast, data_predict[(inst + 1):(inst + horizon)], col = "blue", lwd = 2)
abline(v = monthly_avg$Date[inst], lty = "dashed", col = "gray")
legend("topleft", legend = c("Données", "Historique prédit", "Prévision"),
       col = c("black", "red", "blue"), pch = c(20, NA, NA), lty = c(NA, 1, 1))

#############################################
# 7. Validation : Block Cross-Validation
#############################################
#### Block Cross Validation

# Nombre de blocs
Nblock <- 30
borne_block <- seq(1, nrow(monthly_avg) + 1, length.out = Nblock + 1) %>% floor()
block_list <- list()
l <- length(borne_block)

# Création des blocs
for (i in 2:l) {
  block_list[[i-1]] <- seq(borne_block[i-1], borne_block[i] - 1)
}
block_list <- block_list[lengths(block_list) > 0]  # Retirer les blocs vides

# Fonction pour calculer les résidus pour un bloc donné
block_res <- function(data, alpha, block, smooth.type = "simple") {
  if (smooth.type == "simple") {
    forecast <- expSmooth(data[-block], alpha)
  } else if (smooth.type == "double") {
    forecast <- DoubleExpSmooth(data[-block], alpha)$smooth
  } else {
    stop("Type de lissage non reconnu")
  }
  
  full_forecast <- rep(NA, length(data))
  
  # Insérer les prévisions pour le bloc aux bons indices
  full_forecast[block] <- forecast[seq_along(block)]
  
  # Calculer et retourner les résidus
  return(data[block] - full_forecast[block])
}


# Paramètres
alpha <- 0.5
smooth.type <- "double"

# Calcul des résidus pour tous les blocs
Block_residuals <- lapply(block_list, function(block) {
  block_res(monthly_avg$total_rides, alpha, block, smooth.type)
}) %>% unlist()

# Calcul du MSE moyen
mean_mse <- mean(Block_residuals^2, na.rm = TRUE)
print(paste("Erreur quadratique moyenne : ", mean_mse))