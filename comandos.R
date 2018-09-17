#Header----
# Análises da
# Distribuição espacial de morcegos em gradiente de paisagem na transição com Cerrado
# Parte da tese de doutorado de Nayara Carvalho | Ecologia e Conservação | UFMS

## Preparação
source("funcoes_adicionais.r")
source("pacotes_utilizados.r")
### Para exibir vírgula como separador decimal nos resultados:
options(OutDec = ",")


## Importar os dados usando o pacote readr, salvando como dadosjo
dadosjo <- read_csv("dadosjo.csv", 
                    col_types = cols(DATA = col_date(format = "%B %d, %Y"), 
                                     HORA = col_time(format = "%H%M")))
View(dadosjo) #Os dados estão apresentados em uma data.frame em que cada linha equivale a um morcego
table(dadosjo$ESPECIE) #abundância por espécie | Os números que aparecem como nome das espécies são recapturas.
dim(dadosjo)
nredes <- aggregate(dadosjo$REDE, list(dadosjo$LOCAL, dadosjo$DATA), range)
range(nredes$x[,2]-nredes$x[,1]+1, na.rm = T)

## Matriz de espécies por observação (rede/local/data): abund.rede.local.data a seguir
names(dadosjo) #variáveis no banco de dados (detalhes no arquivo de metadados: INSERIR O NOME DO ARQUIVO)
attach(dadosjo) #para acessar as variáveis em dadosjo mais facilmente
abund.rede.local.data <- aggregate(N, list(LAT, LONG, REDE, LOCAL, DATA), length)
abund.rede.local.data #número de morcegos por observação (rede/local/data) | n = 303
dim(abund.rede.local.data)
colnames(abund.rede.local.data) <- c("LAT","LONG", "REDE", "LOCAL", "DATA", "nindiv")
colnames(abund.rede.local.data) 
observ <- 1:nrow(abund.rede.local.data)
observ #nomeando (ou numerando) as 303 observações (rede/local/data)
rowSums(table(abund.rede.local.data$DATA, observ)) #número de observações por data
rowSums(table(abund.rede.local.data$LOCAL, observ)) #número de observações por local

View(abund.rede.local.data)
abund.observ <- cbind(observ, abund.rede.local.data)
abund.observ #juntando observ com os dados de abundância por observação

### Agora sim a matriz de espécies por observação (spp)
#### Primeiro incluímos as observações e suas respectivas abundâncias na planilha de dados
#### por morcegos (dadosjo)
dados <- merge(abund.observ, dadosjo, all.x = T) #nesse passo 145 morcegos são excluídos porque 
                                                 #informações sobre coodenadas, rede, local ou data
                                                 #estão faltando (NA em dadosjo)
dados <- dados[order(dados$observ), ] #colocando na ordem das observações
head(dados)
names(dados) #variáveis nos dados
dim(dados) #914 morcegos e 69 variáveis

#### Agora usamos essa nova planilha para gerar a matriz observação X espécies (spp)
spp <- table(dados$observ, dados$ESPECIE)
as.matrix(spp)
spp <- as.matrix(spp)[, -(1:3)] #exclui as espécies designadas com números (recapturas)
spp
dim(spp)

##### Aqui juntamos as espécies que precisam de confirmação (indicadas com cf.) com sua respectiva
##### espécie
colnames(spp) #O que fazer com os cf?
spp[, "A. fimbriatus"] <- spp[, "A. cf. fimbriatus"] + spp[, "A. fimbriatus"]
spp[, "A. obscurus"] <- spp[, "A. cf. obscurus"] + spp[, "A. obscurus"]
spp[, "A. planirostris"] <- spp[, "A. cf. planirostris"] + spp[, "A. planirostris"]
spp[, "D. gnoma"] <- spp[, "Dermanura cf. gnoma"] + spp[, "D. gnoma"]
spp[, "P. lineatus"] <- spp[, "P. cf. lineatus"] + spp[, "P. lineatus"]
spp[, "S. lilium"] <- spp[, "S. cf. lilium"] + spp[, "S. lilium"]
spp[, "S. tildae"] <- spp[, "S. cf. tildae"] + spp[, "S. tildae"]
colnames(spp)
spp <- spp[, -c(2:4, 16, 27, 32:33)]
colnames(spp) #Nenhum Carolia tildae ou M. rufus foi identificado com segurança, por isso permanecem com cf.
which(rowSums(spp)==0) #observação vazia
which(colSums(spp)==0) #nenhuma espécie ausente
spp.c <- spp[-which(rowSums(spp)==0), ] #exclui a observação vazia
dim(spp.c) 
spp.mat <- matrix(spp.c, 302, 29) #vamos precisar dos dados na forma de matriz (302 obs. X 29 spp)
colnames(spp.mat) <- colnames(spp.c)
colnames(spp.mat) #29 espécies em 302 observações (rede / local / data)
View(spp.mat)

## Riqueza total em espécies (curva de acumulação por rarefação)
raref.tot <- specaccum(spp.mat, method = "rarefaction") 
raref.tot #910 morcegos de 29 espécies

### Plotar e salvar as curvas de acumulação de espécies:
jpeg("acumula.jpg", width=20, height=15, units="cm", res=300)
plot(raref.tot, ci.type="polygon", ci.col="gray80", ci.lty=0, 
     xlab="Número de morcegos", 
     ylab="Número de espécies", bty="n", xvar="individuals", xlim = c(0, 1000))
dev.off()
#### Figura. Acumulação de 29 espécies em amostra de 910 morcegos na transição do Cerrado com Pantanal e com Floresta Atlântica. A área em cinza indica o intervalo de confiança de 95% para o número médio de espécies estimado por rarefação.

## Amostragem
View(abund.observ)
amostrag <- abund.observ[-which(rowSums(spp)==0), -7]
dim(amostrag)
colnames(amostrag) #Variáveis que identificam as 302 unidades amostrais (rede / local / data)

## Paisagem
paisag <- aggregate(dados[, 27:34], list(dados$observ), mean)
paisag <- paisag[-which(rowSums(spp)==0), -1]
head(paisag) #Proporção de cada componente da paisagem em escala regional (35 km) e local (5 km)
dim(paisag)
rowSums(paisag) #~200 para cada observação (100 para escala local de 5km e 100 para regional de 35km)

## Paisagem por local e data (105 observações | local / data)
paisag.coord.local.data <-  aggregate(paisag, list(amostrag$LAT, amostrag$LONG, amostrag$LOCAL, amostrag$DATA), mean)
colnames(paisag.coord.local.data)[1:4] <- c("LAT", "LONG", "LOCAL", "DATA")
head(paisag.coord.local.data)
paisag.l.d <- paisag.coord.local.data[, -(1:4)]
head(paisag.l.d)
dim(paisag.l.d)

## Paisagem por local (74 locais amostrados)
paisag.l <- aggregate(paisag.l.d, list(paisag.coord.local.data$LOCAL), mean) [, -1]
paisag.l #componentes da paisagem em ordem alfabética dos locais

## Espécies por local e data (abundância média)
spp.coord.local.data <- aggregate(spp.mat, list(amostrag$LAT, amostrag$LONG, amostrag$LOCAL, amostrag$DATA), mean)
colnames(spp.coord.local.data)[1:4] <- c("LAT", "LONG", "LOCAL", "DATA")
head(spp.coord.local.data)
dim(spp.coord.local.data)

### Esforço por data de coleta
range(spp.coord.local.data$DATA)
range(unclass(spp.coord.local.data$DATA)) #tempo desde 01/01/1970
hist(spp.coord.local.data$DATA, breaks=36, freq = T) #esforço amostral durante o período de estudo

#### se quiser apresentar como figura jpg
jpeg("esforcoMensal.jpg", width = 15, height = 15, units = "cm", res = 300)
hist(spp.coord.local.data$DATA, breaks=as.Date(seq(16392, 17432, length.out = 36), origin = "1970-01-01"), 
     freq = T, 
     xlab = "Período de amostragem mensal (nov/14 - set/17)",
     ylab = "Número de pontos amostrados", main = NULL)
dev.off()
### Figura. Esforço mensal de coleta de morcegos com redes de neblina na transição do Cerrado com Pantanal e Floresta Atlântica.

## Esforço por local
table(spp.coord.local.data$LOCAL) #Datas por local

## Só as espécies por local e data
spp.l.d <- spp.coord.local.data[, -(1:4)]
head(spp.l.d)
dim(spp.l.d)
which(rowSums(spp.l.d)==0) #nenhuma observação vazia

#### Número de espécies por observação
nspp.l.d <- rowSums(spp.l.d > 0)
summary(nspp.l.d)
hist(nspp.l.d, breaks = 10)
abline(v = 3, lwd = 3, lty = 2) #mediana do número de espécies por observação
#### Em 50% das observações registramos três espécies ou menos.

### Espécies por local (74 locais amostrados)
spp.l <- aggregate(spp.l.d, list(spp.coord.local.data$LOCAL), mean) [, -1]
spp.l #espécies em ordem alfabética dos locais
dim(spp.l)

#### Número de espécies por local
nspp.l <- rowSums(spp.l > 0)
summary(nspp.l)
hist(nspp.l, breaks = 10)
abline(v = 3, lwd = 3, lty = 2) #mediana do número de espécies por observação
#### Em 50% dos locais registramos três espécies ou menos.

# Área de estudo no Google Maps.
local.data <- spp.coord.local.data[, c("LOCAL", "DATA")]
local.data
coord.l.d <- spp.coord.local.data[, c("LONG", "LAT")]
coord.l.d
coord.l <- aggregate(coord.l.d, list(local.data$LOCAL), mean)
coord.l #coordenadas em ordem alfabética dos locais

## Obtendo a imagem com a função get_maps() do pacote ggmap
imagem <- get_map(location = colMeans(coord.l[,-1]), 
                  zoom = 6, 
                 maptype = "satellite")
### Salvando essa imagem
sat.image <- ggmap(imagem) #

### Adicionando os pontos de amostragem (LOCAL)
image.local <- sat.image +
  geom_point(aes(x = LONG, y = LAT), data = coord.l, 
             size = 1, colour = "orange") + 
  ylab("Latitude") + 
  xlab("Longitude")
image.local

### A mesma imagem mas com tamanhos dos pontos proporcional ao número de espécies
image.local <- sat.image +
  geom_point(aes(x = LONG, y = LAT), data = coord.l, 
             size = nspp.l[order(coord.l[, 1])]*.8, colour = "orange", shape = 21) + 
  ylab("Latitude") + 
  xlab("Longitude") +
  annotate("text", x = -53, y = -24, label = "Rio Paraná", color = "white", 
           size = 3)
image.local
ggsave("riqueza.jpg", width=15, height=15, unit="cm", dpi=300)
#### Figura. Distribuição do número de espécies (entre uma e 10) em 105 observações de 74 locais na transição do Cerrado com Pantanall e Floresta Atlântica. O tamanho dos pontos é diretamente proporcional ao número médio de espécies de morcegos capturados por rede de neblina.

#### Pode juntar essa figura com a curva de acumulação...

# Composição de espécies (diversidade Beta) por observação (local e data)
## Primeiro transformamos para abundância relativa
spp.rel <- decostand(spp.l.d, "total")
rowSums(spp.rel) #só pra conferir: todos 1
colMeans(spp.rel) #abundância relativa média (por rede) de cada espécie
colSums(spp.rel>0) #frequência de ocorrência por observação
spp.rel2 <- spp.rel[, which(colSums(spp.rel>0)>=4)] #excluímos espécies com abundância menor que 4 morcegos

## Ordenar os pontos amostrados em cada data usando NMDS pelas diferenças Bray-Curtis em composição de espécies: 
nmds.2d <- metaMDS(spp.rel2, noshare = T, trymax = 1000) #rodamos até conseguir soluções convergentes
nmds.2d

## Shepard plot
stressplot(nmds.2d)
### Diagrama de Shepard para as distâncias obtidas na NMDS acima com R2 = 0,96 para o ajuste não-métrico.
### Ou seja, as distâncias entre pares de observações no plano da ordenação recuperaram 96% da variância na matriz de distâncias Bray-Curtis.

## Gráfico da ordenação
plot(nmds.2d) #prévia do gráfico de ordenação: círculos são as observações e + as espécies
nmds.escores <- scores(nmds.2d)
nmds.escores # NMDS scores ("posição dos pontos no plano da ordenação")

spp.loads <- wascores(nmds.escores, spp.rel2)
spp.loads #Loadings das espécies ("posição das espécies no plano da ordenação"). Indica quanto cada espécies contribuiu para a ordenação das amostras.

### Gráfico editado:
plot(nmds.escores, bty = "n",  col = "gray", 
     xlab = "Dimensão NMDS 1", 
     ylab = "Dimensão NMDS 2")
abline(h=0, v=0, lty=3, col="gray")
text(spp.loads, rownames(spp.loads),
     col = "gray20", cex = .8, font = 3)
#### Figura. Ordenação das amostras de morcegos por análise não-métrica de escalas multidimensionais (NMDS) a partir de distâncias Bray-Curtis pela abundância relativa. 
#### A posição do nome indica quanto cada espécie contribuiu para a ordenação das amostras (círculos).

### Salvar a figura da ordenação
jpeg("nmds.jpg", width=15, height=15, units="cm", res=300)
plot(nmds.escores, bty = "n",  col = "gray",
     xlab = "Dimensão NMDS 1", 
     ylab = "Dimensão NMDS 2")
abline(h=0, v=0, lty=3, col="gray")
text(spp.loads*1.12, rownames(spp.loads),
     col = "gray20", cex = .8, font = 3)
dev.off()

## PCA da paisagem por observação (local e data)
### PCA de correlação (argumento scale=TRUE)
land.pca <- rda(paisag.l.d, scale=TRUE)
land.pca 
summary(land.pca) # Os dois primeiros eixos da PCA recuperaram 62% da variância na matriz de paisagem, 
                  # sendo 32% no primeiro eixo e 30% no segundo.

#### escores dos locais
site.pca <- land.pca$CA$u
site.pca
#### escores dos componentes da paisagem (loadings)
paisag.pca <- land.pca$CA$v
paisag.pca
#### Eigenvalues | comprimento do eixos | proporcional a variância
(ev <- land.pca$CA$eig)

### Gráfico da PCA (biplot) 
plot(site.pca[, 1:2]*5.4, bty = "n", #essa correção de 5.4 é o "General scaling constant of scores" em summary(land.pca)
     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), col = "gray30",
     xlab = "PCA 1", ylab = "PCA 2")
text(paisag.pca[1:4, 1] * ev[1], paisag.pca[1:4, 2] * ev[2],
     c("Floresta", "Cultivo", "Solo", "Água"))
arrows(0, 0, paisag.pca[1:4, 1] * ev[1]*.8, 
       paisag.pca[1:4, 2] * ev[2]*.8, length = .1)
text(paisag.pca[5:8, 1] * ev[1], paisag.pca[5:8, 2] * ev[2],
     c("Floresta", "Cultivo", "Solo", "Água"), 
     col = "gray60")
arrows(0, 0, paisag.pca[5:8, 1] * ev[1] * .8, 
       paisag.pca[5:8, 2] * ev[2] * .8, length = .1, 
       col = "gray60")
abline(v = 0, h = 0, lty = 3, lwd = .8, col = "gray")

### Bubbles para cada componente da paisagem
quartz() #em Windows use windows()
par(mfrow = c(1, 2), mar = c(5, 4, 1, .1))
plot(coord.l.d[, 1], coord.l.d[, 2],cex.axis = 0.8, pch = 21, 
     col = "white", bg = "light green", ylim = c(-25, -16), xlim = c(-57, -47),
     cex = 4 * paisag.l.d$CULTIVO_NDVI_35/max(paisag.l.d$CULTIVO_NDVI_35),
     xlab = "Longitude", ylab = "Latitude", bty = "n", 
     main = "Regional (raio 35 km)", cex.main = 1)
points(coord.l.d[, 1]+.3, coord.l.d[, 2]-.3, pch=21, col="white", 
       bg="darkgreen", 
       cex=4 * paisag.l.d$FLORESTADA_NDVI_35 / max(paisag.l.d$FLORESTADA_NDVI_35))
points(coord.l.d[, 1]+.6, coord.l.d[, 2]-.6, pch=21, col="white", 
       bg="blue", 
       cex=4 * paisag.l.d$AGUA_SOLO_NDVI_35/max(paisag.l.d$AGUA_SOLO_NDVI_35))
points(coord.l.d[, 1]+.9, coord.l.d[, 2]-.9, pch=21, col="white", 
       bg="bisque", 
       cex=4 * paisag.l.d$CULTIVO_SOLO_NDVI_35/max(paisag.l.d$CULTIVO_SOLO_NDVI_35))
plot(coord.l.d[, 1], coord.l.d[, 2],cex.axis = 0.8, pch = 21, 
     col = "white", bg = "light green", ylim = c(-25, -16), xlim = c(-57, -47),
     cex = 4 * paisag.l.d$CULTIVO_NDVI_5/max(paisag.l.d$CULTIVO_NDVI_5), yaxt = "n",
     xlab = "Longitude", ylab = "", bty = "n", main = "Local (raio 5 km)", cex.main = 1)
points(coord.l.d[, 1]+.3, coord.l.d[, 2]-.3, pch=21, col="white", 
       bg="darkgreen", 
       cex=4 * paisag.l.d$FLORESTADA_NDVI_5 / max(paisag.l.d$FLORESTADA_NDVI_5))
points(coord.l.d[, 1]+.6, coord.l.d[, 2]-.6, pch=21, col="white", 
       bg="blue", 
       cex=4 * paisag.l.d$AGUA_SOLO_NDVI_5/max(paisag.l.d$AGUA_SOLO_NDVI_5))
points(coord.l.d[, 1]+.9, coord.l.d[, 2]-.9, pch=21, col="white", 
       bg="bisque", 
       cex=4 * paisag.l.d$CULTIVO_SOLO_NDVI_5/max(paisag.l.d$CULTIVO_SOLO_NDVI_5))
legend("bottomright", c("Cultivo", "Floresta", "Água", "Solo"), 
       pch = 19, col = c("light green", "dark green","blue", "bisque"),
       bty = "n", text.col = c("light green", "dark green","blue", "bisque"))

### Tudo junto
jpeg("paisagem_PCA_local_data.jpg", width = 20, height = 22.5, units = "cm", res = 300)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = T), 
       heights = c(4, 6))
par(mar = c(.1, 4, 1, 4))
plot(site.pca[, 1:2]*4.9, bty = "n", 
     col = "gray30", xlim = c(-2, 1.5), ylim = c(-1.5, 1.5),
     xlab = "", ylab = "", xaxt = "n", 
     yaxt = "n", main = "PCA | Paisagem")
text(paisag.pca[1:4, 1] * ev[1], paisag.pca[1:4, 2] * ev[2],
     c("Floresta", "Cultivo", "Solo", "Água"))
arrows(0, 0, paisag.pca[1:4, 1] * ev[1]*.8, 
       paisag.pca[1:4, 2] * ev[2]*.8, length = .1)
text(paisag.pca[5:8, 1] * ev[1], paisag.pca[5:8, 2] * ev[2],
     c("Floresta", "Cultivo", "Solo", "Água"), 
     col = "gray60")
arrows(0, 0, paisag.pca[5:8, 1] * ev[1] * .8, 
       paisag.pca[5:8, 2] * ev[2] * .8, length = .1, 
       col = "gray60")
abline(v = 0, h = 0, lty = 3, lwd = .8, col = "gray")

par(mar = c(4, 4, 1, 1))
set.seed(25121971)
plot(jitter(coord.l.d[, 1], 2500), coord.l.d[, 2],cex.axis = 0.8, pch = 21, 
     col = "white", bg = "light green", ylim = c(-25, -16), xlim = c(-57, -45),
     cex = 4 * paisag.l.d$CULTIVO_NDVI_35/max(paisag.l.d$CULTIVO_NDVI_35),
     xlab = "Longitude", ylab = "Latitude", bty = "n", 
     main = "Regional (raio 35 km)", cex.main = 1, col.main = "gray60")
set.seed(25121971)
points(jitter(coord.l.d[, 1], 2500)+.3, coord.l.d[, 2]-.3, pch=21, col="white", 
       bg="darkgreen", 
       cex=4 * paisag.l.d$FLORESTADA_NDVI_35 / max(paisag.l.d$FLORESTADA_NDVI_35))
set.seed(25121971)
points(jitter(coord.l.d[, 1], 2500)+.6, coord.l.d[, 2]-.6, pch=21, col="white", 
       bg="blue", 
       cex=4 * paisag.l.d$AGUA_SOLO_NDVI_35/max(paisag.l.d$AGUA_SOLO_NDVI_35))
set.seed(25121971)
points(jitter(coord.l.d[, 1], 2500)+.9, coord.l.d[, 2]-.9, pch=21, col="white", 
       bg="bisque", 
       cex=4 * paisag.l.d$CULTIVO_SOLO_NDVI_35/max(paisag.l.d$CULTIVO_SOLO_NDVI_35))
set.seed(25121971)
plot(jitter(coord.l.d[, 1], 2500), coord.l.d[, 2],cex.axis = 0.8, pch = 21, 
     col = "white", bg = "light green", ylim = c(-25, -16), xlim = c(-57, -45),
     cex = 4 * paisag.l.d$CULTIVO_NDVI_5/max(paisag.l.d$CULTIVO_NDVI_5), yaxt = "n",
     xlab = "Longitude", ylab = "", bty = "n", main = "Local (raio 5 km)", cex.main = 1)
set.seed(25121971)
points(jitter(coord.l.d[, 1], 2500)+.3, coord.l.d[, 2]-.3, pch=21, col="white", 
       bg="darkgreen", 
       cex=4 * paisag.l.d$FLORESTADA_NDVI_5 / max(paisag.l.d$FLORESTADA_NDVI_5))
set.seed(25121971)
points(jitter(coord.l.d[, 1], 2500)+.6, coord.l.d[, 2]-.6, pch=21, col="white", 
       bg="blue", 
       cex=4 * paisag.l.d$AGUA_SOLO_NDVI_5/max(paisag.l.d$AGUA_SOLO_NDVI_5))
set.seed(25121971)
points(jitter(coord.l.d[, 1], 2500)+.9, coord.l.d[, 2]-.9, pch=21, col="white", 
       bg="bisque", 
       cex=4 * paisag.l.d$CULTIVO_SOLO_NDVI_5/max(paisag.l.d$CULTIVO_SOLO_NDVI_5))
legend("bottomright", c("Cultivo", "Floresta", "Água", "Solo"), 
       pch = 19, col = c("light green", "dark green","blue", "bisque"),
       bty = "n", text.col = c("light green", "dark green","blue", "bisque"))
dev.off()
#### Figura. Ordenação de amostras por análise de componentes principais (PCA) para correlação entre componentes da paisagem (definidos como a proporção da área de cultivo, floresta, água ou solo exposto).
#### As amostras foram obtidas em 105 observações, entre nov/2014 e set/2017, de 74 locais em escala local (buffer com 5 km de raio) e regional (35 km). 
#### O tamanho do pontos é diretamente proporcional a área ocupada pelo componente da paisagem ao redor do ponto amostral.



######################################################################
# Decomposição da diversidade [índice de entropia de Rao (Rao 1982)] #
#                       Adonis-permanova                             #
######################################################################
## Detalhes no script RaoRel.R ou RaoAdo.R (de Bello et al. 2011)

### Diversidade observada (alfa, beta e gama)
Raospp.mat <- RaoRel(sample=t(spp.l.d), dfunc=vegdist(decostand(t(spp.l.d), "total")), dphyl=NULL, 
                weight=T, Jost=F, structure=NULL)
Raospp.mat$TD$Alpha #diversidade de Simpson em cada observação#
witRao<-Raospp.mat$TD$Mean_Alpha
witRao #ALFA (diversidade de Simpson)
betRao<-Raospp.mat$TD$Beta_add #diversidade Beta#
betRao #BETA
totRao<-Raospp.mat$TD$Gamma #diversidade Gamma#
totRao #GAMA
(betRao+witRao)==totRao #BETA + ALFA = GAMA

### Diversidade esperada (alfa, beta e gama)
RaoPerm<-RaoAdo(sample=t(spp.l.d), dfunc=vegdist(decostand(t(spp.l.d), "total")), dphyl=NULL, weight=T, Jost=F, structure=NULL)
RaoPerm$TD$Gamma
RaoPerm$TD$Mean_Alpha
RaoPerm$TD$Beta_add

##Teste de permutação: Permanova para saber se a variação em diversidade entre observações
##(local e data) é diferente da esperada ao acaso (função adonis2 do pacote vegan)
adonis(decostand(spp.l.d, "total") ~ local.data$LOCAL:local.data$DATA)

### A diversidade beta é maior do que a esperada ao acaso segundo a Permanova acima.

#Tabela. Partição da diversidade de morcegos na trasição do Cerrado com Pantanal e Florsta Atlântica.
#Fonte de variação      diversidade obs.	%	    diversidade esp.	%	    
#Interna (alfa)          			     0.436  50.6	           0.519  60.2	  obs < esp	
#Entre observações (beta)			     0.425	49.4	           0.343	39.8	  obs > esp	
#Total	(gama)		                 0.863	                 0.862		      obs = esp	

# Enfim as sementes
sementes <- dados[, 47:ncol(dados)]
head(sementes)
dim(sementes) #sementes de 23 espécies de plantas em uma amostra de 914 morcegos
sum(rowSums(sementes)>0) #223 desses defecaram sementes
223/914*100 #24%

## Sementes por morcegos
dim.temp <- table(dados$ESPECIE, sementes[, 1])
dim.temp
dim(dim.temp)
tab.seed <- matrix(NA, nrow = nrow(dim.temp), ncol = ncol(sementes))
for(i in 1:ncol(sementes)){
  tab.seed[, i] <- table(dados$ESPECIE, sementes[, i])[, 2]
}
tab.seed
rownames(tab.seed) <- rownames(dim.temp)
colnames(tab.seed) <- colnames(sementes)
tab.seed

tab.seed2 <- tab.seed[-which(rowSums(tab.seed)==0), ]
rownames(tab.seed2) 
#### juntar os cf. com as respectivas espécies
tab.seed2["A. fimbriatus", ] <- tab.seed2["A. cf. fimbriatus", ] + tab.seed2["A. fimbriatus", ]
tab.seed2["A. planirostris", ] <- tab.seed2["A. cf. planirostris", ] + tab.seed2["A. planirostris", ]
tab.seed2["P. lineatus", ] <- tab.seed2["P. cf. lineatus", ] + tab.seed2["P. lineatus", ]
tab.seed2["S. tildae", ] <- tab.seed2["S. cf. tildae", ] + tab.seed2["S. tildae", ]
tab.seed2 <- tab.seed2[-c(2, 3, 10, 12), ]
tab.seed2

tab.seed3 <- tab.seed2[which(rowSums(tab.seed2)>3), ]
tab.seed3 <- tab.seed3[, which(colSums(tab.seed3)>3)]
rowSums(tab.seed3)
colSums(tab.seed3)
tab.seed3
tab.seed.rel <- matrix(decostand(t(tab.seed3), "total"), 12, 8)
tab.seed.rel
rownames(tab.seed.rel) <- colnames(tab.seed3)
colnames(tab.seed.rel) <- rownames(tab.seed3)
colnames(tab.seed.rel) <- c("Artibeus_fimbriatus", 
                            "Artibeus_lituratus", 
                            "Artibeus_planirostris", 
                            "Carolia_perspicillata",
                            "Glossophaga_soricina", 
                            "Plathyrrinus_lineatus", 
                            "Sturnira_lilium",
                            "Sturnira_tildae")
rownames(tab.seed.rel) <- c("C. pachystachya",
                            "F. insipida",
                            "F. obtusifolia",
                            "M. tinctoria",
                            "Piper",
                            "P. regnellii",
                            "P. aduncum",
                            "P. amalago",
                            "Psidium",
                            "Rubiaceae",
                            "Solanum 1",
                            "Solanum 2")

### Ordenação das espécies de sementes pela frequência nas espécies de morcegos
nmds.seed <- metaMDS(tab.seed.rel, k = 1, trymax = 200, noshare = T)
stressplot(nmds.seed)

jpeg("sem_mor.jpg", width = 20, height = 25, units = "cm", res = 300)
layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(1.5, 8.5))
par(mar=c(4, 4.1, .1, 12.1))
plot(1:12,rep(0,12),type="n",xaxt="n",yaxt="n", bty="n", xlab="", ylab="")
axis(1, at = 1:12, las=2, cex.axis = 1.2, font = 3,
     labels = rownames(tab.seed.rel)[order(scores(nmds.seed))])
poncho(data.frame(tab.seed.rel), gradient = scores(nmds.seed), 
       xlab = "", ylab = "Frequência relativa", show.grad = NULL, 
       cex.lab = 1.5, cex.species = 1, col = "gray")
dev.off()
### Figura. Espécies de sementes ordenadas por análise não-métrica de escalas multidimensionais (NMDS) pela frequência relativa nas fezes de morcegos.

## Espécies de sementes por rede
sementes.rede.local.data <- aggregate(sementes, 
                                      list(dados$LAT, dados$LONG, 
                                           dados$REDE, 
                                           dados$LOCAL, 
                                           dados$DATA), 
                                      sum)
dim(sementes.rede.local.data)
sementes.rede.local.data
### excluir observação ausente (sem nenhum morcego)
sementes.c <- sementes.rede.local.data[-which(rowSums(spp)==0), ]
head(sementes.c)
dim(sementes.c)
sementes.c <- sementes.c[, -(1:5)]
sort(colSums(sementes.c))
sort(colSums(sementes.c))/302

## Espécies de sementes por local e data (frequência média)
semente.coord.local.data <- aggregate(sementes.c, list(amostrag$LAT, amostrag$LONG, amostrag$LOCAL, amostrag$DATA), mean)
colnames(semente.coord.local.data)[1:4] <- c("LAT", "LONG", "LOCAL", "DATA")
head(semente.coord.local.data)

### Só as espécies de sementes
sementes.l.d <- semente.coord.local.data[, -(1:4)]
head(sementes.l.d)
dim(sementes.l.d)
which(rowSums(sementes.l.d)==0) #observações vazias
sum(rowSums(sementes.l.d)==0)/105 #Em 48% das observações não obtivemos sementes
sementes.l.d.c <- sementes.l.d[-which(rowSums(sementes.l.d)==0), ] #exclui observações vazias

### Espécies de semente por local (74 locais amostrados)
sementes.l <- aggregate(sementes.l.d, list(semente.coord.local.data$LOCAL), mean) [, -1]
sementes.l #espécies em ordem alfabética dos locais

#### Número de espécies por local
nsementes.l <- rowSums(sementes.l > 0)
summary(nsementes.l)
hist(nsementes.l, breaks = 9)
abline(v = 1, lwd = 3, lty = 2) #mediana do número de espécies por observação


#### A variação em diversidade de sementes depende das espécies de morcegos?
##### Primeiro excluímos as observações vazias quanto a sementes da planilha de morcegos
spp.l.d.c <- spp.l.d[-which(rowSums(sementes.l.d)==0), ]
dim(spp.l.d.c)
dim(sementes.l.d.c)

##### Depois transformamos para valores relativos (abundância para morcegos e frequência para sementes)
semente.rel <- decostand(sementes.l.d.c, "total")
spp.c.rel <- decostand(spp.l.d.c, "total")

# Análise espacial
## Transformação dos dados
### Hellinger para os morcegos
spp.hel <- decostand(spp.l, "hel")
dim(spp.hel)

### Distância métrica para as coordenadas
coord2 <- geoXY(coord.l[, 3], coord.l[, 2])
coord2 #distâncias euclidianas desde os pontos mais ao sul (X) e mais ao leste (Y)

# Análise de mapas de autovetores de Moran baseados em distância (dbMEM) com a função quickMEM() de Daniel Borcard (Borcard et. al. 2011)
bats.dbmem.quick <- quickMEM(spp.hel, coord2)
## Essa análise não detectou nenhuma estrutura espacial na comunidade de morcegos.

# A comunidade de morcegos responde a paisagem independentemente do local e período de amostragem?
## Primeiro olhamos para a correlação da paisagem com a comunidade
eixos.pca <- site.pca[, 1:2] # Dois primeiros eixos da PCA da paisagem por local e data
e.pca1 <- eixos.pca[, 1] + abs(min(eixos.pca[, 1])) # tranformação para o eixo 1 manter a escala somente com valores positivos 
e.pca2 <- eixos.pca[, 2] + abs(min(eixos.pca[, 2])) # tranformação para o eixo 2 manter a escala somente com valores positivos
pca.loads <- wascores(nmds.escores, cbind(e.pca1, e.pca2), expand = T)
pca.loads # Escores dos componentes da paisagem ponderados pela média para a ordenação das amostras pelas espécies de morcegos (NMDS).

### Salvar a figura da ordenação com a variação na paisagem
jpeg("nmds_paisag.jpg", width=15, height=15, units="cm", res=300)
plot(nmds.escores, bty = "n",  col = "gray",
     xlab = "Dimensão NMDS 1", 
     ylab = "Dimensão NMDS 2")
abline(h=0, v=0, lty=3, col="gray")
text(spp.loads*1.12, rownames(spp.loads),
     col = "gray50", cex = .8, font = 3)
plot(envfit(nmds.escores, eixos.pca), col = "black", cex = .7)
dev.off()

summary(manova(nmds.escores~site.pca[, 1]+site.pca[, 2]+local.data$LOCAL+local.data$DATA))
## A estrutura da comunidade de morcegos respondeu a variação na paisagem, tanto quanto ao 
## padrão definido pela substituição de cobertura florestal por cultivo (PCA 1), quanto ao
## padrão definido pela substituição da proporção de áreas com água por solo exposto (PCA 2).
## A comunidade de morcegos também diferiu entre os locais de amostragem, reforçando a 
## importância da diversidade beta, como evidenciado com a permanova (~50% da variância na 
## diversidade total de morcegos deveu-se a variação entre observações | diversidade beta).

# Agora voltamos as sementes
sementes.l.d.c  #Frequência média de ocorrência de sementes (número de morcegos em que sementes ocorreram nas fezes)
rownames(sementes.l.d.c)
rowSums(sementes.l.d.c)
semente.rel  #Frequência relativa de ocorrência de sementes
rownames(semente.rel)
dim(semente.rel)
rowSums(sementes.l.d.c > 0) #número de espécies de sementes
sementes.l.d.c[which(rowSums(sementes.l.d.c>0)==1),] #observações com uma espécie de semente
colSums(sementes.l.d.c > 0) #frequência de ocorrência das espécies de sementes
seed.util.nmds <- sementes.l.d.c[-which(rowSums(sementes.l.d.c>0)==1),] #observações com mais de uma espécie de semente
dim(seed.util.nmds)
seed.rel.util.nmds <- semente.rel[-which(rowSums(sementes.l.d.c>0)==1),]
dim(seed.rel.util.nmds) #agora com frequência relativa

spp.l.d.c #Abundância média de morcegos (número de morcegos por rede)
dim(spp.l.d.c)
rownames(spp.l.d.c)
spp.l.d.c[which(rowSums(spp.l.d.c>0)==1),]
rownames(spp.c.rel) #Abundância relativa das espécies de morcegos
spp.c.rel.seed <- spp.c.rel[-which(rowSums(sementes.l.d.c>0)==1),]
spp.c.rel.seed$`D. rotundus`
dim(spp.c.rel.seed)
which(colSums(spp.c.rel.seed)>0) #espécies de morcegos em que registramos mais de uma espécie de semente 
nmds.bats.seed <- metaMDS(spp.c.rel.seed, noshare = T)
stressplot(nmds.bats.seed)

paisag.l.d.c <- paisag.l.d[-which(rowSums(sementes.l.d)==0),]
paisag.l.d.c #Porcentagem da área ocupada por cada um dos componentes da paisagem em escala local (5 km) e regional (35 Km)
dim(paisag.l.d.c)
rownames(paisag.l.d.c)
paisag.l.d.c2 <- paisag.l.d.c[-which(rowSums(sementes.l.d.c>0)==1),]
dim(paisag.l.d.c2)
pca.land.seed <- prcomp(paisag.l.d.c2, scale. = T)
summary(pca.land.seed)
eixos.pca.seed <- pca.land.seed$x[, 1:2]


l.d.c <- local.data[-which(rowSums(sementes.l.d)==0),]
l.d.c #Local e data das amostras de sementes
dim(l.d.c)
l.d.c2 <- l.d.c[-which(rowSums(sementes.l.d.c>0)==1),]

### Em 55 observações (local / data) registramos sementes de 23 espécies (ou morfoespécies) nas fezes de
### morcegos de 29 espécies

seed.nmds <- metaMDS(seed.rel.util.nmds, noshare = T)
stressplot(seed.nmds)
seeds.escores <- scores(seed.nmds)
seeds.loads <- wascores(seeds.escores, seed.rel.util.nmds)
seeds.loads

jpeg("nmds_sementes_blue2.jpg", width=18, height=15, units="cm", res=300)
plot(seeds.escores, bty = "n",  col = "gray75",
     xlab = "Dimensão NMDS 1", 
     ylab = "Dimensão NMDS 2")
abline(h=0, v=0, lty=3, col="gray")
spp.bats.seeds <- spp.c.rel.seed[, which(colSums(spp.c.rel.seed>0)>3)]
plot(envfit(seeds.escores, spp.bats.seeds), adj = 0, font=3, arrow.mul=1.5, col = "blue", cex = .7)
seeds.loads["Piper sp2", ] <- c(-0.01407567-.05,  0.141063604+.05)
text(seeds.loads, rownames(seeds.loads), adj = 0.3,
     cex = .7, font = 3)
dev.off()
## Figura. Ordenação de amostras (círculos cinza) por análise não-métrica de escalas multidimensionais (NMDS) pela composição de espécies de sementes nas fezes de morcegos.
## As setas indicam as correlações da abundância relativa das espécies mais frequentes (> 3 ocorrências) de morcegos com o plano de ordenação.

## Permanova 
colSums(seed.rel.util.nmds)
colSums(spp.c.rel.seed)
spp.ad <- spp.c.rel.seed[-which(colSums(spp.c.rel.seed)==0)]
semente.adon <- adonis(as.matrix(seed.rel.util.nmds) ~ 
                       as.matrix(spp.ad) + 
                       as.matrix(eixos.pca.seed), 
                       permutations = 999, set.seed(25121971))
semente.adon #73% da variância nas diferenças quanto a composição de sementes foi explicada pela composição de espécies de morcegos independentemente da paisagem
