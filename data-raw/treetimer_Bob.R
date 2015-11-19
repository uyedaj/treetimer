require(geiger)
require(rotl)
require(foreach)
require(phytools)
require(aRbor)
require(doParallel)
require(phyndr)
require(devtools)
require(roxygen2)
require(distory)
rm(list=ls(all=TRUE))
load_all()

getTTOL(list.only=TRUE)
ttol <- getTTOL(2)
#save(ttol, file= "./data/2.TTOL_all_smoothed.Rda")
#ottTable <- getOttIds(ttol$tip.label, ncores=12)
#save(ottTable, file="./data/ottTable_2.all_smoothed.Rda")
load("./data/ottTable_2.all_smoothed.Rda")
ottTable[,1] <- unname(sapply(gsub(" ", "_", ottTable[,1]), simpleCap))
td <- make.treedata(ttol, ottTable)
ttolObject <- filter(td, !is.na(ott_id))
tree <- ttolObject$phy
dat <- ttolObject$dat
TH <- max(branching.times(tree))
times <- c(TH-1, 1000, 750, 400, 300, 250, 200, 150, 100, 75, 50, 25, 10, 5, 1)
rm(td, ttol, ottTable)
#lineages <- extractLineages(dat$ott_id, ncores=8)
#save(lineages, file="./data/ottLineages.Rda")
load("./data/ottLineages.Rda")
ttolObject$lineages <- lineages
taxalist <- c("Drosophila_willistoni", "Felis_silvestris", "Homo_sapiens", "Aedes_vexans",
              "Gallus_gallus", "Tinca_tinca", "Cotinis_nitida","Turdus xanthorhynchus",
              "Phrynobatrachus_dispar", "Plethodon_shermani","Driloleirus_macelfreshi", "Driloleirus_americanus", "Bulimulus_darwini")
bobslist <- c("Ellychnia_corrusca", "Cicindela_campestris", "Photinus_pyralis", "Photinus_marginellus",
              "Melaleuca_quinquenervia", "Melaleuca_decora", "Pyractonema_angulata", "Chauliognathus_pennsylvanicus",
              "Chauliognathus", "Tabanus_catenatus", "Grammoptera", "Lactuca_sativa",
              "Photuris_lucicrescens", "Plecia", "Bombus", "Tabanus_abactor", "Gopherus_agassizii",
              "Gopherus polyphemus", "Alligator_mississippiensis", "Capra aegagrus", "Meleagris_gallopavo",
              "Columba", "Homo_sapiens", "Bos_taurus", "Panthera_leo", "Gallus_gallus", "Cercopithecus_aethiops",
              "Tinca_tinca", "Canis_familiaris", "Mus_musculus", "Phoca_vitulina", "Sturnus_vulgaris", "Felis_domestica",
              "Panorpa_communis", "Apis_mellifera", "Atrichopogon_geminus", "Atrichopogon_levis", "Cantharis_bilineatus", "Chrysops",
              "Neoaliturus_haematoceps", "Cotinis_nitida", "Aedes_sollicitans", "Diabrotica_undecimpunctata",
              "Culex_annulus", "Eriocheir_sinensis", "Liriodendron_tulipifera", "Tabanus_gladiator", "Bidens",
              "Ixodes_pacificus", "Dalbulus_maidis", "Photuris_pennsylvanicus", "Leptinotarsa_decemlineata",
              "Leucoma_salicis", "Tabanus_lineola", "Tabanus_nigrovittatus", "Haemaphysalis_leporispalustris", "Monobia_quadridens",
              "Hybomitra_opaca", "Penaeus_vannamei", "Macrosteles_fascifrons", "Pachydiplax_longipennis", "Drosophila_willistoni",
              "Aedes_sticticus", "Aedes_vexans", "Trigonotylus_ruficornis", "Cocos_nucifera", "Prunus", "Eristalis_arbustorum",
              "Tabanus_abdominalis", "Culex_tritaeniorhynchus", "Haematopota")

bobOtts <- getOttIds(bobslist, ncores=12)

ttPhynd <- phyndrTTOL(ttolObject, bobslist, times, prune=TRUE, ncores=5)
newTimetree <- phyndr_sample(ttPhynd$phyndr)
plot(newTimetree)
length(newTimetree$tip.label)/length(taxalist)*100
synth_tree <- tol_induced_subtree(ttPhynd$otts$ott_id)
synth_tree$tip.label <- strip_ott_ids(synth_tree$tip.label)
synth_tree <- compute.brlen(synth_tree, power = 2)
plot(synth_tree)
conTree <- congruify.phylo(newTimetree, synth_tree, scale="PATHd8")
edgecols <- rep(1, nrow(conTree$phy$edge))
edgecols[distinct.edges(conTree$phy, newTimetree)] <- 2
plot(conTree$phy, edge.color=edgecols)

save(list(ttPhynd, conTree), file="./data/bobsTreePhynd.Rda")
