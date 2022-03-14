## produce a random, non-ultrametric tree
library(ape)
library(RRphylo)
library(visdat)

# Exemplo package RRPhylo
rtree(100)->tree
rtree
tree

####
?setBM
setBM(tree,nY=1,type="brown")-> y

x <- setBM(tree,nY=1,type="drift")

z <- setBM(tree,nY=1,type="trend")

#Observe que setBM difere de fastBM em que os fenótipos produzidos são verificados quanto à existência de uma tendência temporal no fenótipo. O usuário pode especificar se deseja dados sem tendência (opção "marrom"), fenótipos com tendência no tempo (opção "drift") ou fenótipos cuja variância aumenta/diminui exponencialmente ao longo do tempo, consistente com a existência de uma tendência na taxa de evolução (opção "tendência"). Neste último caso, o usuário pode indicar a intensidade da tendência (aplicando diferentes valores de es), e se ela deve ocorrer após uma determinada proporção da altura da árvore (portanto, um determinado ponto de volta no tempo, especificado pelo argumento t .mudança). Árvores em setBM são tratadas como não ultramétricas. Se uma árvore ultramétrica é alimentada à função, setBM altera ligeiramente os comprimentos das folhas multiplicando aleatoriamente metade das folhas por 1 * 10e-3, a fim de torná-la não ultramétrica.

## run phylogenetic ridge regression
RRphylo(tree,y)->RR
hist(RR$rates)

#$tree, a versão totalmente dicotômica do argumento tree; $tip.path, a matriz L; $node.path, a matriz L′; $taxas, o vetor das taxas evolutivas; $ace, o vetor de estimativas do estado ancestral; $y.estimates, o vetor de estados de ponta estimados; e $lambda, que é o fator de penalização ajustado

#O cálculo do fator de penalização usa a função optL (embutida no RRphylo), para encontrar a estimativa de máxima verossimilhança de λ minimizando a variação da taxa da raiz da árvore em direção às pontas

# Grandes valores absolutos de λ são consistentes com sinal filogenético muito alto (o que significa que espécies filogeneticamente próximas tenderão a ter fenótipos muito semelhantes, veja a vinheta para explicação completa). Descobrimos empiricamente que os valores de λ de 0 a 1 são consistentes com o modelo de evolução do movimento browniano



# setBM produz um vetor ou matriz fenotípica com uma tendência desejada no fenótipo ao longo do tempo ou na taxa evolutiva ao longo do tempo
search.shift(RR, status.type="clade", filename=getwd(), "SSauto") 

search.shift(RR, status.type = "sparse",node = NULL, state= NULL, cov = NULL, nrep = 1000, f = NULL,filename=paste(tempdir()),"test_phy")

#auto.recognize="yes",,covariate= "FALSE"

# Como o RRphylo atribui uma taxa específica a cada ramo da árvore, é viável aplicar indiferentemente quando os diferentes regimes de taxa pertencem a clados distintos ou a várias espécies não relacionadas ao longo da filogenia. Nós nos referimos a essas duas situações distintas como condições distribuídas em nível de “clado” e “esparsas” (filogeneticamente). Sob a condição “esparsa”, uma reconstrução em nível de clado é teoricamente possível, dada a estimativa apropriada dos estados. No entanto, preferimos evitar inferências sobre estados ancestrais, uma vez que os métodos atualmente disponíveis baseiam-se na suposição frequentemente violada de que o traço de interesse evolui a uma taxa uniforme em toda a árvore (King & Lee, 2015). A função R que fornecemos para localizar mudanças de taxa evolutiva e testar hipóteses específicas sobre a taxa anexada a diferentes clados ou estados de ponta é search.shift. Ele pega um objeto produzido pelo RRphylo e pode ser usado para localizar automaticamente os deslocamentos, testar diferentes clados onde se presume que regimes de taxa distintos se apliquem ou testar diferenças de taxa entre diferentes estados de ponta

## Not run: 
data("DataOrnithodirans")
DataOrnithodirans$treedino->treedino
DataOrnithodirans$massdino->massdino
DataOrnithodirans$statedino->statedino
View(DataOrnithodirans)

View(massdino)
RRphylo(tree=treedino,y=massdino)->dinoRates
class(treedino)
class(anura_phy)
amphibia_phy <- ape::multi2di(anura_phy)
# Case 1. Without accounting for the effect of a covariate

# Case 1.1 "clade" condition
# with auto-recognize
search.shift(RR=dinoRates,status.type="clade",
             filename=paste(tempdir(), "SSauto", sep="/"))

# testing two hypothetical clades
search.shift(RR=dinoRates,status.type="clade",node=c(696,746),
             filename=paste(tempdir(), "SSclade", sep="/"))

# Case 1.2 "sparse" condition
# testing the sparse condition.
search.shift(RR=dinoRates,status.type= "sparse",state=statedino,
             filename=paste(tempdir(), "SSsparse", sep="/"))


# Case 2. Accounting for the effect of a covariate

# Case 2.1 "clade" condition
search.shift(RR=dinoRates,status.type= "clade",cov=massdino,
             filename=paste(tempdir(), "SSclade_cov", sep="/"))

# Case 2.2 "sparse" condition
search.shift(RR=dinoRates,status.type="sparse",state=statedino,cov=massdino,
             filename=paste(tempdir(), "SSstate_cov", sep="/"))

## End(Not run)

###Eventualmente, search.shift pode testar ainda mais o efeito de uma covariável nos valores de taxa e usar os resíduos da regressão de taxa versus covariável (em vez dos valores de taxa absoluta ajustados pelo RRphylo) para contrastar taxas entre diferentes ramificações. É importante notar que a covariável deve ter o mesmo comprimento que o vetor de taxa (ou seja, duas vezes o número de pontas, menos uma). Se a covariável for os próprios fenótipos para os quais as taxas são calculadas, o usuário deve apenas especificar o fenótipo como o argumento “cov”. Caso contrário, se um fenótipo diferente for selecionado como covariável, as taxas para este segundo fenótipo devem ser submetidas primeiro ao RRphylo e, em seguida, as estimativas ancestrais e os valores de ponta coletados nesta ordem e alimentados para search.shift como o objeto “cov”

# Thu Mar 10 16:14:06 2022 ------------------------------
# Testando com os dados de anfíbios
library(phytools)
library(picante)

# Pegando filogenia anfíbios
setwd("E:/OneDrive/2019 - Moroti/Tese/Capitulo 1/RMarkdown/phylogeny")
anura.phy <- read.tree("phy_amphibia.tre")

# Pegando traits anfíbios
setwd("E:/OneDrive/2019 - Moroti/Tese/Capitulo 1/RMarkdown/trait_data")
anura.trait <- read.csv2("traits_amphibia.csv")
# selecionando colunas
amphibia_full_traits <- anura.trait[ , 1:2]

# nomenclatura
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp, "Boana_claresignata", "Bokermannohyla_claresignata")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp, "Boana_clepsydra", "Bokermannohyla_clepsydra")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Proceratophrys_salvatori", "Odontophrynus_salvatori")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Scinax_v_signatus", "Scinax_v-signatus")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Scinax_x_signatus", "Scinax_x-signatus")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Boana_bandeirantes", "Hypsiboas_bandeirantes")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Boana_poaju", "Hypsiboas_poaju")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Boana_paranaiba", "Hypsiboas_paranaiba")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Physalaemus_nattereri", "Eupemphix_nattereri")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Agalychnis_aspera", "Hylomantis_aspera")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Agalychnis_granulosa", "Hylomantis_granulosa")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Lithobates_palmipes", "Rana_palmipes")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Pithecopus_azureus", "Phyllomedusa_azurea")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Pithecopus_ayeaye", "Phyllomedusa_ayeaye")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Pithecopus_nordestinus", "Phyllomedusa_nordestina")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Pithecopus_rohdei", "Phyllomedusa_rohdei")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Pithecopus_rohdei", "Phyllomedusa_rohdei")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Aparasphenodon_ararapa", "Aparasphenodon_arapapa")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Ololygon_tripui", "Scinax_tripui")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Ololygon_cosenzai", "Scinax_cosenzai")
amphibia_full_traits$sp <- str_replace_all(amphibia_full_traits$sp,"Ololygon_strigilata", "Scinax_strigilatus")
amphibia_full_traits$sp<-str_replace_all(amphibia_full_traits$sp,"Pithecopus_hypochondrialis", "Phyllomedusa_hypochondrialis")
amphibia_full_traits$sp<-str_replace_all(amphibia_full_traits$sp,"Pithecopus_megacephalus", "Phyllomedusa_megacephala")

# Removendo missing data
short.amphibia.traits <- remove_missing(short.amphibia.traits, vars=names(short.amphibia.traits))
names(short.amphibia.traits) 

# passando nome para as colunas para poder comparar
test <- xtabs( ~ Boh

###
# Precisamos remover nos traits
rem.col.phy <- c("Brachycephalus_sulfuratus", "Crossodactylus_fransciscanus","Crossodactylus_timbuhy","Dendropsophus_bromeliaceus","Eleutherodactylus_bilineatus","Melanophryniscus_milanoi","Melanophryniscus_xanthostomus","Proceratophrys_mantiqueira","Proceratophrys_pombali","Pseudopaludicola_atragula","Pseudopaludicola_pocoto","Scinax_strigilatus","Trachycephalus_typhonius","Chiasmocleis_alagoana","Crossodactylus_werneri","Ololygon_tupinamba","Phyllodytes_megatympanum","Physalaemus_atim","Proceratophrys_phyllostoma","Proceratophrys_tupinamba","Pseudopaludicola_jaredi","Rhinella_sebbeni")

traits.amphibia.phy <- short.amphibia.traits[!short.amphibia.traits$sp %in% rem.col.phy,]

# Agora precisamos cortar a filogenia com as sp nos traits
anura_phy <- prune.sample(test, anura.phy)

# Agora precisamos transformar os dados para entrar na função RRphylo
# primeiro colocar os nomes das espécies nas colunas e o tamanho nas linhas
#traits_amphibia_phy <- traits.amphibia.phy %>%
#  group_by(sp) %>%
#  mutate(Grp = row_number()) %>%
#  pivot_wider(names_from = "sp", values_from = "Body_size_mm") %>%
#  select(-Grp) 
# traits_amphibia_phy <- traits_amphibia_phy[-2,]

# Precisamos converter esse objeto para númerico
svl <- log10(as.numeric(traits.amphibia.phy$Body_size_mm))
names(svl) <- traits.amphibia.phy$sp

# Taxa de evolução
rates.anura <- RRphylo(tree= amphibia_phy, y=svl)
rates_anura <- rates.anura$rates[546:1091,]
