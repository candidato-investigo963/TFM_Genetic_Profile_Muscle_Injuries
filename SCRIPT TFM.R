###############################################################################
# SCRIPT TFM : Análisis Curvas ROC y Regresión Logística
###############################################################################

# AUTOR: Pablo Rafael Pombero Hurtado

# TFM: Perfil genético de rendimiento muscular y su relación con la incidencia 
# y gravedad de lesiones musculoesqueléticas en futbolistas no profesionales

# AVISO: Este script se ha utilizado con la variable binaria FIBRASmusc, 
# pero puede usarse también con cualquier otra variable binaria.

### Librerías utilizadas 

library(dplyr)
library(readxl)
library(writexl)
library(ROCR)
library(pROC)
library(ggplot2)
library(tidyr)


### Cálculo del promedio del TGS


# Leer base de datos
datos <- read_excel("C:/Users/Win10/OneDrive/Escritorio/MASTER/PRACTICAS/Lesiones_Genética.xlsx")

# Eliminar datos ausentes (NA)
datos_limpios <- datos %>% filter(!is.na(`TGS rendimiento`) & !is.na(FIBRASmusc))

# Agrupar por la variable de estudio y calcular media y desviación estándar del TGS
resumen <- datos_limpios %>%
  group_by(FIBRASmusc) %>%
  summarise(
    media_TGS = mean(`TGS rendimiento`),
    sd_TGS = sd(`TGS rendimiento`)
  ) %>%
  mutate(Grupo = ifelse(FIBRASmusc == 0, "No lesionados", "Lesionados"))

# Distribución del TGS según los estados de la variable binaria
tapply(datos$`TGS rendimiento`, datos$FIBRASmusc, mean, na.rm = TRUE)

# Gráfico de barras
ggplot(resumen, aes(x = Grupo, y = media_TGS, fill = Grupo)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = media_TGS - sd_TGS, ymax = media_TGS + sd_TGS),
                width = 0.2, size = 1) +
  scale_fill_manual(values = c("black", "gray")) +
  labs(
    title = "Puntuación Genética Promedio (TGS) según lesión de fibras musculares",
    y = "TGS (a.u.)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Aplicar t-student para ver diferencias significativas en las medias
t.test(`TGS rendimiento` ~ FIBRASmusc , data = datos)




### Análisis mediante CURVA ROC


# Filtrar solo no profesionales
datos <- subset(datos, Profesional == 0)

# Eliminar datos ausentes (NA)
datos <- na.omit(datos[, c("SEXO", "TGS rendimiento", "FIBRASmusc")])

# Filtrar por sexo
hombres <- subset(datos, SEXO == 0)
mujeres <- subset(datos, SEXO == 1)

# ===== HOMBRES ROC =====
roc_h <- roc(hombres$FIBRASmusc, hombres$`TGS rendimiento`, ci = TRUE)
auc_h <- as.numeric(auc(roc_h))
ci_h <- ci.auc(roc_h)
n1_h <- sum(hombres$FIBRASmusc == 1)
n0_h <- sum(hombres$FIBRASmusc == 0)
Q1_h <- auc_h / (2 - auc_h)
Q2_h <- 2 * auc_h^2 / (1 + auc_h)
se_auc_h <- sqrt((auc_h*(1 - auc_h) + (n1_h - 1)*(Q1_h - auc_h^2) + (n0_h - 1)*(Q2_h - auc_h^2)) / (n1_h * n0_h))
z_h <- (auc_h - 0.5) / se_auc_h
p_h <- round(2 * (1 - pnorm(abs(z_h))), 4)
coords_h <- coords(roc_h, x = "best", best.method = "youden", ret = c("threshold", "specificity", "sensitivity"))
thresh_h <- round(as.numeric(coords_h["threshold"]), 2)
spec_h <- as.numeric(coords_h["specificity"])
sens_h <- as.numeric(coords_h["sensitivity"])

# ===== MUJERES ROC =====
roc_m <- roc(mujeres$FIBRASmusc, mujeres$`TGS rendimiento`, ci = TRUE)
auc_m <- as.numeric(auc(roc_m))
ci_m <- ci.auc(roc_m)
n1_m <- sum(mujeres$FIBRASmusc == 1)
n0_m <- sum(mujeres$FIBRASmusc == 0)
Q1_m <- auc_m / (2 - auc_m)
Q2_m <- 2 * auc_m^2 / (1 + auc_m)
se_auc_m <- sqrt((auc_m*(1 - auc_m) + (n1_m - 1)*(Q1_m - auc_m^2) + (n0_m - 1)*(Q2_m - auc_m^2)) / (n1_m * n0_m))
z_m <- (auc_m - 0.5) / se_auc_m
p_m <- round(2 * (1 - pnorm(abs(z_m))), 4)
coords_m <- coords(roc_m, x = "best", best.method = "youden", ret = c("threshold", "specificity", "sensitivity"))
thresh_m <- round(as.numeric(coords_m["threshold"]), 2)
spec_m <- as.numeric(coords_m["specificity"])
sens_m <- as.numeric(coords_m["sensitivity"])



# ===== GRÁFICO ROC HOMBRES =====
roc_df_h <- data.frame(
  fpr = as.numeric(1 - roc_h$specificities),
  tpr = as.numeric(roc_h$sensitivities)
)

ggplot(roc_df_h, aes(x = fpr, y = tpr)) +
  geom_line(color = "red", size = 1.3) +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1.1) +
  labs(
    title = "Curva ROC – TGS rendimiento y lesión muscular (hombres)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0,
           label = paste0("AUC ", round(auc_h, 3), 
                          " (", round(ci_h[1], 3), " - ", round(ci_h[3], 3), 
                          ") | p = ", p_h),
           size = 5.5, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "#f2f2f2", color = NA),
    panel.background = element_rect(fill = "#f2f2f2", color = NA),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )


ggplot(roc_df_h, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  annotate("point", x = 1 - spec_h, y = sens_h, color = "red", size = 3) +
  annotate("text", x = 1 - spec_h + 0.05, y = sens_h, 
           label = paste0("Umbral: ", thresh_h), hjust = 0, size = 5, color = "red") +
  labs(title = "Curva ROC - Masculino",
       subtitle = paste0("AUC = ", round(auc_h, 3), 
                         " | IC 95%: ", round(ci_h[1], 3), "-", round(ci_h[3], 3),
                         " | p = ", p_h),
       x = "1 - Especificidad", y = "Sensibilidad") +
  theme_minimal(base_size = 14)


# ===== GRÁFICO ROC MUJERES =====
roc_df_m <- data.frame(
  fpr = as.numeric(1 - roc_m$specificities),
  tpr = as.numeric(roc_m$sensitivities)
)

ggplot(roc_df_m, aes(x = fpr, y = tpr)) +
  geom_line(color = "red", size = 1.3) +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1.1) +
  labs(
    title = "Curva ROC – TGS rendimiento y lesión muscular (mujeres)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  annotate("text", x = 0.95, y = 0.05, hjust = 1, vjust = 0,
           label = paste0("AUC ", round(auc_m, 3), 
                          " (", round(ci_m[1], 3), " - ", round(ci_m[3], 3), 
                          ") | p = ", p_m),
           size = 5.5, fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "#f2f2f2", color = NA),
    panel.background = element_rect(fill = "#f2f2f2", color = NA),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )


ggplot(roc_df_m, aes(x = fpr, y = tpr)) +
  geom_line(color = "darkorange", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  annotate("point", x = 1 - spec_m, y = sens_m, color = "red", size = 3) +
  annotate("text", x = 1 - spec_m + 0.05, y = sens_m, 
           label = paste0("Umbral: ", thresh_m), hjust = 0, size = 5, color = "red") +
  labs(title = "Curva ROC - Femenino",
       subtitle = paste0("AUC = ", round(auc_m, 3), 
                         " | IC 95%: ", round(ci_m[1], 3), "-", round(ci_m[3], 3),
                         " | p = ", p_m),
       x = "1 - Especificidad", y = "Sensibilidad") +
  theme_minimal(base_size = 14)



### Regresion logística binaria

# Filtrar solo no profesionales
datos <- subset(datos, Profesional == 0)

# Eliminar datos ausentes (NA)
datos <- na.omit(datos[, c("SEXO", "TGS rendimiento", "FIBRASmusc")])


# ===== HOMBRES =====
hombres <- subset(datos, SEXO == 0 )
hombres$grupo <- ifelse(hombres$`TGS rendimiento` > 66.23, 1, 0) # Usar el umbral de corte obtenido en la curva ROC
modelo_h <- glm(FIBRASmusc ~ grupo, data = hombres, family = binomial)
OR_h <- round(exp(coef(modelo_h)[2]), 3)
IC_h <- round(exp(confint(modelo_h)[2, ]), 3)
pval_h <- round(summary(modelo_h)$coefficients[2, 4], 4)

# ===== MUJERES =====
mujeres <- subset(datos, SEXO == 1)
mujeres$grupo <- ifelse(mujeres$`TGS rendimiento` > 66.23, 1, 0)
modelo_m <- glm(FIBRASmusc ~ grupo, data = mujeres, family = binomial)
OR_m <- round(exp(coef(modelo_m)[2]), 3)
IC_m <- round(exp(confint(modelo_m)[2, ]), 3)
pval_m <- round(summary(modelo_m)$coefficients[2, 4], 4)


# Resultados
cat("===== HOMBRES =====\n")
cat("OR:", OR_h, "| IC 95%:", IC_h[1], "-", IC_h[2], "| p =", pval_h, "\n\n")
cat("===== MUJERES =====\n")
cat("OR:", OR_m, "| IC 95%:", IC_m[1], "-", IC_m[2], "| p =", pval_m, "\n")



