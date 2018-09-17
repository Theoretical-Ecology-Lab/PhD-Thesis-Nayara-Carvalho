# Análises da
# Distribuição espacial de morcegos em gradiente de paisagem na transição com Cerrado
# Parte da tese de doutorado de Nayara Carvalho | Ecologia e Conservação | UFMS

## Objetos
dadosjo <- read_csv("dadosjo.csv", col_types = cols(DATA = col_date(format = "%B %d, %Y"), HORA = col_time(format = "%H%M")))
abund.rede.local.data <- aggregate(dadosjo$N, list(dadosjo$LAT, dadosjo$LONG, dadosjo$REDE, dadosjo$LOCAL, dadosjo$DATA), length)
colnames(abund.rede.local.data) <- c("LAT","LONG", "REDE", "LOCAL", "DATA", "nindiv")
observ <- 1:nrow(abund.rede.local.data)
abund.observ <- cbind(observ, abund.rede.local.data)
dados <- merge(abund.observ, dadosjo, all.x = T)
dados <- dados[order(dados$observ), ] #colocando na ordem das observações
spp <- table(dados$observ, dados$ESPECIE)
spp <- as.matrix(spp)[, -(1:3)] #exclui as espécies designadas com números (recapturas)
spp[, "A. fimbriatus"] <- spp[, "A. cf. fimbriatus"] + spp[, "A. fimbriatus"]
spp[, "A. obscurus"] <- spp[, "A. cf. obscurus"] + spp[, "A. obscurus"]
spp[, "A. planirostris"] <- spp[, "A. cf. planirostris"] + spp[, "A. planirostris"]
spp[, "D. gnoma"] <- spp[, "Dermanura cf. gnoma"] + spp[, "D. gnoma"]
spp[, "P. lineatus"] <- spp[, "P. cf. lineatus"] + spp[, "P. lineatus"]
spp[, "S. lilium"] <- spp[, "S. cf. lilium"] + spp[, "S. lilium"]
spp[, "S. tildae"] <- spp[, "S. cf. tildae"] + spp[, "S. tildae"]
spp <- spp[, -c(2:4, 16, 27, 32:33)]
spp.c <- spp[-which(rowSums(spp)==0), ]
spp.mat <- matrix(spp.c, 302, 29) #vamos precisar dos dados na forma de matriz (302 obs. X 29 spp)
colnames(spp.mat) <- colnames(spp.c)
raref.tot <- specaccum(spp.mat, method = "rarefaction") 
amostrag <- abund.observ[-which(rowSums(spp)==0), -7]
paisag <- aggregate(dados[, 27:34], list(dados$observ), mean)
paisag <- paisag[-which(rowSums(spp)==0), -1]
paisag.coord.local.data <-  aggregate(paisag, list(amostrag$LAT, amostrag$LONG, amostrag$LOCAL, amostrag$DATA), mean)
colnames(paisag.coord.local.data)[1:4] <- c("LAT", "LONG", "LOCAL", "DATA")
paisag.l.d <- paisag.coord.local.data[, -(1:4)]
paisag.l <- aggregate(paisag.l.d, list(paisag.coord.local.data$LOCAL), mean) [, -1]
spp.coord.local.data <- aggregate(spp.mat, list(amostrag$LAT, amostrag$LONG, amostrag$LOCAL, amostrag$DATA), mean)
colnames(spp.coord.local.data)[1:4] <- c("LAT", "LONG", "LOCAL", "DATA")
spp.l.d <- spp.coord.local.data[, -(1:4)]
nspp.l.d <- rowSums(spp.l.d > 0)
spp.l <- aggregate(spp.l.d, list(spp.coord.local.data$LOCAL), mean) [, -1]
nspp.l <- rowSums(spp.l > 0)
local.data <- spp.coord.local.data[, c("LOCAL", "DATA")]
coord.l.d <- spp.coord.local.data[, c("LONG", "LAT")]
coord.l <- aggregate(coord.l.d, list(local.data$LOCAL), mean)
imagem <- get_map(location = colMeans(coord.l[,-1]), zoom = 6, maptype = "satellite")
sat.image <- ggmap(imagem) 
image.local <- sat.image +
  geom_point(aes(x = LONG, y = LAT), data = coord.l, size = nspp.l[order(coord.l[, 1])]*.8, colour = "orange", shape = 21) + 
  ylab("Latitude") + 
  xlab("Longitude") +
  annotate("text", x = -53, y = -24, label = "Rio Paraná", color = "white", size = 3)
spp.rel <- decostand(spp.l.d, "total")
spp.rel2 <- spp.rel[, which(colSums(spp.rel>0)>=4)] #excluímos espécies com abundância menor que 4 morcegos
nmds.2d <- metaMDS(spp.rel2, noshare = T, trymax = 200) #rodamos até conseguir soluções convergentes
nmds.escores <- scores(nmds.2d)
spp.loads <- wascores(nmds.escores, spp.rel2)
land.pca <- rda(paisag.l.d, scale=TRUE)
site.pca <- land.pca$CA$u
paisag.pca <- land.pca$CA$v
ev <- land.pca$CA$eig
Raospp.mat <- RaoRel(sample=t(spp.l.d), dfunc=vegdist(decostand(t(spp.l.d), "total")), dphyl=NULL, weight=T, Jost=F, structure=NULL)
witRao<-Raospp.mat$TD$Mean_Alpha
betRao<-Raospp.mat$TD$Beta_add #diversidade Beta#
totRao<-Raospp.mat$TD$Gamma #diversidade Gamma#
RaoPerm<-RaoAdo(sample=t(spp.l.d), dfunc=vegdist(decostand(t(spp.l.d), "total")), dphyl=NULL, weight=T, Jost=F, structure=NULL)
sementes <- dados[, 47:ncol(dados)]
sementes.rede.local.data <- aggregate(sementes, list(dados$LAT, dados$LONG, dados$REDE, dados$LOCAL, dados$DATA), sum)
sementes.c <- sementes.rede.local.data[-which(rowSums(spp)==0), ]
sementes.c <- sementes.c[, -(1:5)]
semente.coord.local.data <- aggregate(sementes.c, list(amostrag$LAT, amostrag$LONG, amostrag$LOCAL, amostrag$DATA), mean)
colnames(semente.coord.local.data)[1:4] <- c("LAT", "LONG", "LOCAL", "DATA")
sementes.l.d <- semente.coord.local.data[, -(1:4)]
sementes.l.d.c <- sementes.l.d[-which(rowSums(sementes.l.d)==0), ] #exclui observações vazias
sementes.l <- aggregate(sementes.l.d, list(semente.coord.local.data$LOCAL), mean) [, -1]
nsementes.l <- rowSums(sementes.l > 0)
spp.l.d.c <- spp.l.d[-which(rowSums(sementes.l.d)==0), ]
semente.rel <- decostand(sementes.l.d.c, "total")
spp.c.rel <- decostand(spp.l.d.c, "total")
spp.hel <- decostand(spp.l, "hel")
coord2 <- geoXY(coord.l[, 3], coord.l[, 2])
eixos.pca <- site.pca[, 1:2] # Dois primeiros eixos da PCA da paisagem por local e data
e.pca1 <- eixos.pca[, 1] + abs(min(eixos.pca[, 1])) # tranformação para o eixo 1 manter a escala somente com valores positivos 
e.pca2 <- eixos.pca[, 2] + abs(min(eixos.pca[, 2])) # tranformação para o eixo 2 manter a escala somente com valores positivos
pca.loads <- wascores(nmds.escores, cbind(e.pca1, e.pca2), expand = T)
sementes.l.d.c  #Frequência média de ocorrência de sementes (número de morcegos em que sementes ocorreram nas fezes)
seed.util.nmds <- sementes.l.d.c[-which(rowSums(sementes.l.d.c>0)==1),] #observações com mais de uma espécie de semente
seed.rel.util.nmds <- semente.rel[-which(rowSums(sementes.l.d.c>0)==1),]
spp.c.rel.seed <- spp.c.rel[-which(rowSums(sementes.l.d.c>0)==1),]
nmds.bats.seed <- metaMDS(spp.c.rel.seed, noshare = T)
paisag.l.d.c <- paisag.l.d[-which(rowSums(sementes.l.d)==0),]
paisag.l.d.c2 <- paisag.l.d.c[-which(rowSums(sementes.l.d.c>0)==1),]
pca.land.seed <- prcomp(paisag.l.d.c2, scale. = T)
eixos.pca.seed <- pca.land.seed$x[, 1:2]
l.d.c <- local.data[-which(rowSums(sementes.l.d)==0),]
l.d.c2 <- l.d.c[-which(rowSums(sementes.l.d.c>0)==1),]
seed.nmds <- metaMDS(seed.rel.util.nmds, noshare = T)
seeds.escores <- scores(seed.nmds)
seeds.loads <- wascores(seeds.escores, seed.rel.util.nmds)
spp.ad <- spp.c.rel.seed[-which(colSums(spp.c.rel.seed)==0)]
semente.adon <- adonis(as.matrix(seed.rel.util.nmds) ~ as.matrix(spp.ad) + as.matrix(eixos.pca.seed), permutations = 999, set.seed(25121971))

