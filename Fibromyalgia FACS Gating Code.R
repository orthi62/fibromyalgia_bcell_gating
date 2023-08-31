#Install packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowCore")

require(tidyverse)
require(flowCore)
require(flowClust)
require(ggcyto)
require(gridExtra)

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowAI")

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("flowGate")

library(flowAI)
library(flowGate)
library(cytolib)
library(cytoinstaller)
library(flowStats)
library(flowDensity)
library(flowViz)
library(flowWorkspaceData)


#List files and parameters from flowset.

load("C:/BMD1.RData")

sampleNames(BMD1_flowset)
colnames(BMD1_flowset) 

#Obtaining markers associated to each dye:

key <- c("$P7S", "$P8S", "$P9S", "$P10S", "$P11S", "$P12S", "$P13S", "$P14S", "$P15S", "$P16S", "$P17S", "$P18S", "$P19S", "$P20S", "$P21S", "$P22S", "$P23S", "$P24S", "$P25S")
keyword(BMD1_flowset, keyword = key)

#or

markernames(BMD1_flowset)

#Check Excel Sheet


#COMPENSATION USING INTRINSIC SPILL MATRIX

##COMPENSATE BY PRECALCULATED MATRIX

BMD1_spill_pre <- spillover(BMD_patientcomb_compensated_trans[[46]])
BMD1_spill_pre
BMD1_spillmatrix_pre <- BMD1_spill_pre[[1]]

BMD1_compensated_pre <- compensate(BMD1_flowset, BMD1_spillmatrix_pre)


#ALTERNATIVE: R-CALCULATION OF SPILL MATRIX

#Create a flowSet of Compensation Control frames

comp_control <- flowSet_to_list(BMD1_flowset[18:33])
comp_controlfs <- as(comp_control, "flowSet")
class(comp_controlfs)

colnames(comp_controlfs)
sampleNames(comp_controlfs)

pData(comp_controlfs)


#name the flowframes in the list to order

names(comp_control) <- c("APC-A", "APC-Cy7-A", "BB515-A", "BUV395-A", "BUV661-A", "BUV737-A", "BV421-A", "BV605-A", "BV650-A", "BV711-A", "BV786-A", "PE-A", "PE-CF594-A", "PE-Cy7-A", "PerCP-A", "UNSTAINED")

class(comp_control) 
#returns a list

#convert list to flowset

comp_controlfs <- as(comp_control, "flowSet")
comp_controlfs

pData(comp_controlfs)

colnames(comp_controlfs)

#At this point there are more channels than there are samples (Alexa 700, BV510 and PE-Cy5 are empty channels as there are no compensation controls for these).

#Remove these channels while keeping the FSC-A and SSC-A channels for computing the spillover matrix.

comp_controlfs2 <- comp_controlfs[,c(1, 4, 7:9, 11:12, 14:22, 24)]
colnames(comp_controlfs2)
compBMD1_spill <- spillover(comp_controlfs2, unstained = "UNSTAINED", patt = "-A", fsc = "FSC-A", ssc = "SSC-A", stain_match = "regexpr")

#stain_match = "regexpr" matches the names to the right channels based on patterns.

compBMD1_spill

#COMPENSATION

#Create flowSet without Compensation Controls

BMD1_nocontrol <- flowSet_to_list(BMD1_flowset[c(1:17, 34:66)])
BMD1_nocontrolfs <- as(BMD1_nocontrol, "flowSet")

class(BMD1_nocontrolfs)

#Compensate flowSet

BMD1_compensated <- compensate(BMD1_nocontrolfs, compBMD1_spill)
BMD1_compensated


#CLEANUP

#must have flowAI and call library(flowAI)

BMD1_compensated_clean <- flow_auto_qc(BMD1_compensated)
BMD1_compensated_clean

#Check if the cleaned up plots look alright.

autoplot(BMD1_compensated_clean[[30]])
ggcyto(BMD1_compensated_clean[[30]], aes(x="FSC-A", y = "SSC-A")) + geom_hex(bins=256)+ggcyto_par_set(limits = list(x = c(0, 1e+05), y= c(0, 1e+05)))


#TRANSFORMATION

logicle <- function(x){
  estimateLogicle(x, colnames(x)[7:24])
}

BMD1_trans_matrix <- fsApply(BMD1_compensated_clean, logicle)

BMD1_compensated_trans <- transform (BMD1_compensated_clean, BMD1_trans_matrix)


#GATING

BMD3_auto_gs <- GatingSet(BMD3_compensated_clean)
class(BMD3_auto_gs)

ggcyto(BMD1_compensated_trans[[46]], aes(x = "FSC-A", y = "SSC-A")) + geom_hex(bins = 256)





#1) TIME GATING:

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs)

ggcyto(BMD1_compensated_clean[[46]], aes(x = "Time")) + geom_density(color = "royalblue3")
g.time <- rectangleGate(filterId = "Time", "Time" = c(4000, 34000))
gs_pop_add(BMD1_auto_gs, g.time, parent = "root", name = "Time")
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[46]], x = "Time", "Time", bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)


#Find count statistics.
BMD1_time_count <- flowCore::filter(BMD1_compensated_trans, g.time)
BMD1_time_count

summary(BMD1_time_count)

ggcyto(BMD1_auto_gs[[46]],aes(x="Time", y= "SSC-A"), subset = "Time")+ geom_hex(bins = 256)
time_plot <- ggcyto(BMD1_auto_gs,aes(x="Time", y= "SSC-A"), subset = "root")+ geom_hex(bins = 2000)+ geom_gate(g.time) + ggcyto_par_set(limits = list(x = c(0,4e+04), y = c(0, 2e+05)))
ggsave(time_plot, file = "/Users/orthi/Documents/time.pdf", height = 20, width = 20)

gs_pop_get_stats(BMD1_auto_gs)

#2) LOOSE GATING:

BMD3_fs_gate_data <- gs_pop_get_data(BMD3_auto_gs)

ggcyto(BMD3_auto_gs[[1]], aes(x = "FSC-A", y = "SSC-A"), subset = "root") + geom_hex(bins = 256)

g.nondebris <- rectangleGate(filterId = "Non Debris","FSC-A"=c(0.25e5,Inf))
gs_pop_add(BMD3_auto_gs, g.nondebris, parent = "root", name = "Non Debris" )
recompute(BMD3_auto_gs)

autoplot(BMD3_auto_gs[[1]], x = "FSC-A", y = "SSC-A", "Non Debris", bins = 256)

gs_get_pop_paths(BMD3_auto_gs)

BMD1_nondebris_count <- flowCore::filter(BMD1_compensated_trans, g.nondebris)
BMD1_nondebris_count

summary(BMD1_nondebris_count)

ggcyto(BMD1_auto_gs[[1]],aes(x="FSC-A", y= "SSC-A"), subset = "Non Debris")+ geom_hex(bins = 256)

gs_pop_get_stats(BMD3_auto_gs)



#3) SINGLET GATING:

BMD3_fs_gate_data <- gs_pop_get_data(BMD3_auto_gs, "Non Debris")

ggcyto(BMD3_auto_gs[[1]], aes(x = "FSC-A", y = "FSC-H"), subset = "Non Debris") + geom_hex(bins = 256)

g.singlet <- polygonGate(filterId = "Singlets", "FSC-A" = c(0.5e4,25e4,25e4,0.5e4),"FSC-H"=c(0e4,18e4,20e4,2e4))
gs_pop_add(BMD3_auto_gs, g.singlet, parent = "Non Debris", name = "Singlets")
recompute(BMD3_auto_gs)

autoplot(BMD3_auto_gs[[1]], x = "FSC-A", y = "FSC-H", "Singlets", bins = 256)

gs_get_pop_paths(BMD3_auto_gs)


#Find count statistics.

BMD1_singlets_count <- flowCore::filter(BMD1_compensated_trans, g.singlet)
BMD1_singlets_count

summary(BMD1_singlets_count)

gs_pop_get_stats(BMD3_auto_gs)

ggcyto(BMD1_auto_gs[[1]],aes(x="FSC-A",y="FSC-H"),subset="root")+geom_hex(bins = 256)+geom_gate("Singlets")+geom_stats(adjust = 0.8)

p5 <- ggcyto(BMD1_auto_gs,aes(x="FSC-A",y="FSC-H"),subset="root")+geom_hex(bins = 1000)+geom_gate("Singlets")+geom_stats(adjust = 0.8)+ggcyto_par_set(limits = list(x = c(0,3.5e5), y = c(0, 3.5e5)))

ggsave(p5, file = "/Users/orthi/Documents/Singlets.pdf", width = 20, height = 20)

ggcyto(BMD1_auto_gs[[46]],aes(x="FSC-A", y= "FSC-H"),subset="Singlets")+ geom_hex(bins = 256)+geom_gate(g.singlet)

ggcyto(BMD1_auto_gs[[46]],aes(x="FSC-A", y= "FSC-H"),subset="Singlets")+ geom_hex(bins = 256)


#4) CELLS OF INTEREST (CD3-/CD56-)

BMD3_fs_gate_data <- gs_pop_get_data(BMD3_auto_gs, "Singlets")

ggcyto(BMD3_auto_gs[[1]], aes(x = "BV605-A", y = "APC-Cy7-A"), subset = "Singlets") + geom_hex(bins = 256)

#To identify minimum and maximum gate parameters for x and y-axis:

CD56_FMO <- ggcyto(BMD1_compensated_trans$FMO_CD56.fcs, aes(x = "BV605-A", y = "SSC-A"), subset = "root") + geom_hex(bins = 256)
ggsave(CD56_FMO, file = "/Users/orthi/Documents/CD56_FMO.pdf", width = 20, height = 20)


gate.CD3_minus_CD56_minus <- quadGate(filterId = "CD3-/CD56-", "BV605-A"=2.2, "APC-Cy7-A"=2)
gs_pop_add(BMD3_auto_gs, gate.CD3_minus_CD56_minus, parent = "Singlets", name = c("CD3+CD56-", "CD3+CD56+", "CD3-CD56+", "CD3-CD56-"))
recompute(BMD3_auto_gs)

autoplot(BMD3_auto_gs[[1]], x = "BV605-A", y = "APC-Cy7-A", c("CD3+CD56-", "CD3+CD56+", "CD3-CD56+", "CD3-CD56-"), bins = 256)

gs_get_pop_paths(BMD3_auto_gs)

#Find count statistics.

BMD1_cells_count <- flowCore::filter(BMD1_compensated_trans, gate.CD3_minus_CD56_minus)
BMD1_cells_count

summary(BMD1_cells_count)

gs_pop_get_stats(BMD2_auto_gs)

ggcyto(BMD1_auto_gs[[46]],aes(x="BV605-A",y="APC-Cy7-A"),subset="Singlets")+geom_hex(bins = 256)+geom_gate(gate.CD3_minus_CD56_minus)+geom_stats(adjust = 0.8)
DUMP_remove <- ggcyto(BMD1_auto_gs,aes(x="BV605-A",y="APC-Cy7-A"), subset="Singlets")+geom_hex(bins = 256)+geom_gate(gate.CD3_minus_CD56_minus)+geom_stats(adjust = 0.8)
ggsave(BMD1_auto_gs, file = "/Users/orthi/Documents/DUMP_remove.pdf", width = 20, height = 20)

ggcyto(BMD1_auto_gs[[46]],aes(x="BV605-A", y= "APC-Cy7-A"),subset="CD3-CD56-")+ geom_hex(bins = 256)+geom_gate(gate.CD3_minus_CD56_minus)

ggcyto(BMD1_auto_gs[[46]],aes(x="BV605-A", y= "APC-Cy7-A"),subset="CD3-CD56-")+ geom_hex(bins = 256)

#5) CD19 GATE

BMD3_fs_gate_data <- gs_pop_get_data(BMD3_auto_gs, "CD3-CD56-")

ggcyto(BMD3_auto_gs[[1]], aes(x = "BB515-A", y = "SSC-A"), subset = "CD3-CD56-") + geom_hex(bins = 256)
gate.CD19 <- rectangleGate(filterId = "CD19+", "BB515-A" = c(2, 4.5), "SSC-A" = c(0, 120000))
gs_pop_add(BMD3_auto_gs, gate.CD19, parent = "CD3-CD56-", name = "CD19+")
recompute(BMD3_auto_gs)

gs_get_pop_paths(BMD3_auto_gs)


#if we need to remove a gate:

gs_pop_remove(BMD1_auto_gs, "/Time/Non Debris/Singlets/CD3-CD56-/CD19-")

#Find count statistics.

BMD1_CD19_count <- flowCore::filter(BMD1_compensated_trans, gate.CD19)
BMD1_CD19_count

summary(BMD1_CD19_count)

gs_pop_get_stats(BMD3_auto_gs)

ggcyto(BMD2_auto_gs[[1]],aes(x="BB515-A",y="SSC-A"),subset="CD3-CD56-")+geom_hex(bins = 256)+geom_gate(gate.CD19)+geom_stats(adjust = 0.8)
ggcyto(BMD1_auto_gs,aes(x="BB515-A",y="SSC-A"),subset="CD3-CD56-")+geom_hex(bins = 256)+geom_gate(gate.CD19)+geom_stats(adjust = 0.8)

ggcyto(BMD1_auto_gs[[46]],aes(x="BB515-A", y= "SSC-A"),subset="CD19+")+ geom_hex(bins = 256)+geom_gate(gate.CD19)

ggcyto(BMD1_auto_gs[[46]],aes(x="BB515-A", y= "SSC-A"),subset="CD19+")+ geom_hex(bins = 256)

p3 = ggcyto(BMD1_auto_gs,aes(x="BB515-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 256)+geom_gate(gate.CD19)+geom_stats(adjust = 0.8)

ggsave(p3, file = "/Users/orthi/Documents/B_cell_percentage_rel_Singlets.pdf", width = 20, height = 20)

#6) CD27/IgD gate

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs, "CD19+")

ggcyto(BMD1_auto_gs[[1]], aes(x = "BV786-A", y = "PE-CF594-A"), subset = "CD19+") + geom_hex(bins = 256)

gate.CD27_IgD <- quadGate(filterId = "CD27/IgD", "BV786-A"=2.2, "PE-CF594-A"=1.0)
gs_pop_add(BMD1_auto_gs, gate.CD27_IgD, parent = "CD19+", name = c("IgD+CD27-", "IgD+CD27+", "IgD-CD27+", "IgD-CD27-"))
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[50]], x = "BV786-A", y = "PE-CF594-A", c("IgD+CD27-", "IgD+CD27+", "IgD-CD27+", "IgD-CD27-"), bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)

#Find count statistics.

BMD1_Bcelltype_count <- filter(BMD1_compensated_trans, gate.CD27_IgD)
BMD1_Bcelltype_count

summary(BMD1_Bcelltype_count)

gs_pop_get_stats(BMD1_auto_gs)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV786-A",y="PE-CF594-A"),subset="root")+geom_hex(bins = 256)+geom_gate(gate.CD27_IgD)+geom_stats(adjust = 0.8)
ggcyto(BMD1_auto_gs,aes(x="BV786-A",y="PE-CF594-A"),subset="CD19+")+geom_hex(bins = 256)+geom_gate(gate.CD27_IgD)+geom_stats(adjust = 0.8)


ggcyto(BMD1_auto_gs[[50]],aes(x="BV786-A", y= "PE-CF594-A"),subset="CD19+")+ geom_hex(bins = 256)+geom_gate(gate.CD27_IgD)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV786-A", y= "PE-CF594-A"),subset="IgD-CD27+")+ geom_hex(bins = 256)


#8) CD38+/- plasmablasts from IgD-CD27+

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs, "IgD-CD27+")

ggcyto(BMD1_auto_gs[[1]], aes(x = "BV421-A", y = "SSC-A"), subset = "IgD-CD27+") + geom_hex(bins = 256)

gate.CD38_plus <- rectangleGate(filterId = "CD38+", "BV421-A" = c(2, 4.6), "SSC-A" = c(20000, 100000))
gs_pop_add(BMD1_auto_gs, gate.CD38_plus, parent = "IgD-CD27+", name = "CD38+")
recompute(BMD1_auto_gs)

gate.CD38_minus <- rectangleGate(filterId = "CD38-", "BV421-A" = c(0.5, 2), "SSC-A" = c(20000, 100000))
gs_pop_add(BMD1_auto_gs, gate.CD38_minus, parent = "IgD-CD27+", name = "CD38-")
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[3]], x = "BV421-A", y = "SSC-A", c("CD38+", "CD38-"), bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)

#Find count statistics.

BMD1_pb_count_Ig_minus_CD27_plus <- filter(BMD1_compensated_trans, gate.CD38_plus)
BMD1_pb_count_Ig_minus_CD27_plus

summary(BMD1_pb_count_Ig_minus_CD27_plus)

BMD1_pbminus_count_Ig_minus_CD27_plus <- filter(BMD1_compensated_trans, gate.CD38_minus)
BMD1_pbminus_count_Ig_minus_CD27_plus

summary(BMD1_pbminus_count_Ig_minus_CD27_plus)

gs_pop_get_stats(BMD1_auto_gs)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV421-A", y="SSC-A"),subset= "IgD-CD27+")+geom_hex(bins = 256)+geom_gate(gate.CD38_plus)+geom_gate(gate.CD38_minus)+geom_stats(adjust = 0.8)
CD38_gates <- ggcyto(BMD1_auto_gs,aes(x="BV421-A",y="SSC-A"),subset="IgD-CD27+")+geom_hex(bins = 256)+geom_gate(gate.CD38_plus)+geom_gate(gate.CD38_minus)+geom_stats(adjust = 0.8)

ggsave(CD38_gates, file = "/Users/orthi/Documents/plasmablast_percentage_rel_memory.pdf", height = 20, width = 20)


#9) CD38+/- plasmablasts Ig+CD27+ 

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs, "IgD+CD27+")

ggcyto(BMD1_auto_gs[[3]], aes(x = "BV421-A", y = "SSC-A"), subset = "IgD+CD27+") + geom_hex(bins = 256)

gs_pop_add(BMD1_auto_gs, gate.CD38_plus, parent = "IgD+CD27+", name = "CD38+p")
recompute(BMD1_auto_gs)

gs_pop_add(BMD1_auto_gs, gate.CD38_minus, parent = "IgD+CD27+", name = "CD38-p")
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[50]], x = "BV421-A", y = "SSC-A", c("CD38+p", "CD38-p"), bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)

#Find count statistics.

BMD1_pbplus_count_IgD_plus_CD27_plus <- filter(BMD1_compensated_trans, gate.CD38_plus)
BMD1_pbplus_count_IgD_plus_CD27_plus

summary(BMD1_pbplus_count_IgD_plus_CD27_plus)

BMD1_pbminus_count_IgD_plus_CD27_plus <- filter(BMD1_compensated_trans, gate.CD38_minus)
BMD1_pbminus_count_IgD_plus_CD27_plus

summary(BMD1_pbminus_count_IgD_plus_CD27_plus)

gs_pop_get_stats(BMD1_auto_gs)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV421-A", y="SSC-A"),subset= "IgD+CD27+")+geom_hex(bins = 256)+geom_gate(gate.CD38_plus)+geom_gate(gate.CD38_minus)+geom_stats(adjust = 0.8)
CD38_gates_IgDp_CD27p <- ggcyto(BMD1_auto_gs,aes(x="BV421-A",y="SSC-A"),subset="IgD+CD27+")+geom_hex(bins = 256)+geom_gate(gate.CD38_plus)+geom_gate(gate.CD38_minus)+geom_stats(adjust = 0.8)

ggsave(CD38_gates_IgDp_CD27p, file = "/Users/orthi/Documents/plasmablast_percentage_rel_memoryIgM.pdf", height = 20, width = 20)


#10) CD10 gating on Ig+CD27-

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs, "IgD+CD27-")

ggcyto(BMD1_auto_gs[[3]], aes(x = "BV711-A", y = "SSC-A"), subset = "IgD+CD27-") + geom_hex(bins = 256)

gate.CD10_plus <- rectangleGate(filterId = "CD10+", "BV711-A" = c(2.2, 3), "SSC-A" = c(0, 100000))
gs_pop_add(BMD1_auto_gs, gate.CD10_plus, parent = "IgD+CD27-", name = "CD10+")
recompute(BMD1_auto_gs)

gate.CD10_minus <- rectangleGate(filterId = "CD10-", "BV711-A" = c(0, 2.2), "SSC-A" = c(0, 100000))
gs_pop_add(BMD1_auto_gs, gate.CD10_minus, parent = "IgD+CD27-", name = "CD10-")
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[50]], x = "BV711-A", y = "SSC-A", c("CD10+", "CD10-"), bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)

#Find count statistics.

BMD1_transitional_count_IgpCD27n_plus <- filter(BMD1_compensated_trans, gate.CD10_plus)
BMD1_transitional_count_IgpCD27n_plus

summary(BMD1_transitional_count_IgpCD27n_plus)

BMD1_transitional_count_IgpCD27n_minus <- filter(BMD1_compensated_trans, gate.CD10_minus)
BMD1_transitional_count_IgpCD27n_minus

summary(BMD1_transitional_count_IgpCD27n_minus)

gs_pop_get_stats(BMD1_auto_gs)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV711-A", y="SSC-A"), subset= "IgD+CD27-")+ geom_hex(bins = 256)+ geom_gate(gate.CD10_plus) + geom_gate(gate.CD10_minus)+geom_stats(adjust = 0.8)
CD10_gates <- ggcyto(BMD1_auto_gs,aes(x="BV711-A",y="SSC-A"),subset="IgD+CD27-")+geom_hex(bins = 256)+geom_gate(gate.CD10_plus)+geom_gate(gate.CD10_minus)+geom_stats(adjust = 0.8)

ggsave(CD10_gates, file = "/Users/orthi/Documents/transitional_percentage_rel_naive.pdf", height = 20, width = 20)


#11) Gate for CD21+ from CD10-

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs, "CD10-")

ggcyto(BMD1_auto_gs[[2]], aes(x = "BUV737-A", y = "SSC-A"), subset = "CD10-") + geom_hex(bins = 256)

gate.CD21_plus <- rectangleGate(filterId = "CD21+", "BUV737-A" = c(1.7, 3.5), "SSC-A" = c(0, 100000))
gs_pop_add(BMD1_auto_gs, gate.CD21_plus, parent = "CD10-", name = "CD21+")
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[50]], x = "BUV737-A", y = "SSC-A", "CD21+", bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)

#Find count statistics.

BMD1_latetrans_anergic_count <- filter(BMD1_compensated_trans, gate.CD21_plus)
BMD1_latetrans_anergic_count

summary(BMD1_latetrans_anergic_count)

gs_pop_get_stats(BMD1_auto_gs)

ggcyto(BMD1_auto_gs[[50]],aes(x="BUV737-A", y="SSC-A"), subset= "CD10-")+ geom_hex(bins = 256)+ geom_gate(gate.CD21_plus) + geom_stats(adjust = 0.8)
CD21_gate <- ggcyto(BMD1_auto_gs,aes(x="BUV737-A",y="SSC-A"),subset="CD10-")+geom_hex(bins = 256)+geom_gate(gate.CD21_plus)+ geom_stats(adjust = 0.8)

ggsave(CD21_gate, file = "/Users/orthi/Documents/anergic_percentage_rel_late_trans_naive.pdf", height = 20, width = 20)

plot(BMD1_auto_gs)

#11) CD38/CD27 gating

BMD1_fs_gate_data <- gs_pop_get_data(BMD1_auto_gs, "IgD-CD27-")

ggcyto(BMD1_auto_gs[[2]], aes(x = "BV786-A", y = "BV421-A"), subset = "IgD-CD27-") + geom_hex(bins = 256)

gate.CD27_CD38 <- rectangleGate(filterId = "CD27-", "BV786-A" = c(0, 2.2), "BV421-A" = c(2, 4))
gs_pop_add(BMD1_auto_gs, gate.CD27_CD38, parent = "IgD-CD27-", name = "CD38+CD27-")
recompute(BMD1_auto_gs)

autoplot(BMD1_auto_gs[[50]], x = "BV786-A", y = "BV421-A", "CD38+CD27-", bins = 256)

gs_get_pop_paths(BMD1_fs_gate_data)
gs_get_pop_paths(BMD1_auto_gs)

#Find count statistics.

BMD1_activatedAtbc_count <- filter(BMD1_compensated_trans, gate.CD27_CD38)
BMD1_activatedAtbc_count

summary(BMD1_activatedAtbc_count)

gs_pop_get_stats(BMD1_auto_gs)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV786-A",y="BV421-A"),subset="IgD-CD27-")+geom_hex(bins = 256)+geom_gate(gate.CD27_CD38)+geom_stats(adjust = 0.8)
CD38_CD27_ <- ggcyto(BMD1_auto_gs,aes(x="BV786-A",y="BV421-A"),subset="IgD-CD27-")+geom_hex(bins = 256)+geom_gate(gate.CD27_CD38)+geom_stats(adjust = 0.8)
CD38C
ggsave(CD38_CD27_, file = "/Users/orthi/Documents/activated_percentage_rel_Atbc.pdf", height = 20, width = 20)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV786-A", y= "BV421-A"),subset="root")+ geom_hex(bins = 256)+geom_gate(gate.CD27_CD38)

ggcyto(BMD1_auto_gs[[50]],aes(x="BV786-A", y= "BV421-A"),subset="CD38+CD27-")+ geom_hex(bins = 256)



#obtain FCS files for each gating population

B3_cell_gate <- cytoset_to_flowSet(gs_pop_get_data(BMD3_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+"))
naive_gate <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27-"))
unswitched_memory <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27+"))
switched_memory <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD-CD27+"))
AtBc <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD-CD27-"))
plasmablast_unswitched <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27+/CD38+p"))
nonpb_unswitched <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27+/CD38-p"))
plasmablast_switched <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD-CD27+/CD38+"))
nonpb_switched <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD-CD27+/CD38-"))
naive_trans <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27-/CD10+"))
naive_nontrans <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27-/CD10-"))
naive_anergic <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD+CD27-/CD10-/CD21+"))
AtBc_activated <- cytoset_to_flowSet(gs_pop_get_data(BMD1_auto_gs, "/Non Debris/Singlets/CD3-CD56-/CD19+/IgD-CD27-/CD38+CD27-"))

write.flowSet(B3_cell_gate, filename = "B_Cell.fcs", outdir = "/Users/orthi/Documents/B_cell_3" )
write.flowSet(naive_gate, filename = "B_Naive.fcs", outdir = "/Users/orthi/Documents/B_Naive")
write.flowSet(unswitched_memory, filename = "B_Unswitched.fcs", outdir = "/Users/orthi/Documents/B_Unswitched" )
write.flowSet(switched_memory, filename = "B_Switched.fcs", outdir = "/Users/orthi/Documents/B_switched" )
write.flowSet(AtBc, filename = "B_Atbc.fcs", outdir = "/Users/orthi/Documents/B_Atbc" )
write.flowSet(plasmablast_switched, filename = "B_pbswitch.fcs", outdir = "/Users/orthi/Documents/B_pbswitch" )
write.flowSet(plasmablast_unswitched, filename = "B_pbunswitch.fcs", outdir = "/Users/orthi/Documents/B_pbunswitch" )
write.flowSet(nonpb_switched, filename = "B_nonpbswitch.fcs", outdir = "/Users/orthi/Documents/B_nonpbswitch" )
write.flowSet(nonpb_unswitched, filename = "B_nonpbunswitch.fcs", outdir = "/Users/orthi/Documents/B_nonpbunswitch")
write.flowSet(naive_trans, filename = "B_Trans.fcs", outdir = "/Users/orthi/Documents/B_Trans" )
write.flowSet(naive_nontrans, filename = "B_NonTrans.fcs", outdir = "/Users/orthi/Documents/B_NonTrans" )
write.flowSet(naive_anergic, filename = "B_Anergic.fcs", outdir = "/Users/orthi/Documents/B_Anergic" )
write.flowSet(AtBc_activated, filename = "B_AtBc_active.fcs", outdir = "/Users/orthi/Documents/B_AtBc_active" )


names <- read.table("/Users/orthi/Documents/B_cell_3/annotation.txt")
names
class(names[1,1])

new_names <- c("Control30_ID_E25_B_cell", "Case27_ID_E29_B_cell", "Case28_ID_E30_B_cell", "Case29_ID_E36_B_cell", "Control31_ID_E38_B_cell", "Control32_ID_E39_B_cell", "Control33_ID_E40_B_cell", "Control34_ID_E44_B_cell", "Control35_ID_E45_B_cell", "Control36_ID_E46_B_cell", "Control37_ID_E47_B_cell", "Control38_ID_E48_B_cell", "Control39_ID_E50_B_cell", "Control40_ID_E52_B_cell", "Control41_ID_E54_B_cell", "Control42_ID_E55_B_cell", "Control42_ID_E59_B_cell", "Case30_ID_E61_B_cell", "Control55_ID_F100_B_cell", "Control56_ID_F110_B_cell", "Case31_ID_F36_B_cell", "Case32_ID_F40_B_cell", "Control44_ID_F43_B_cell", "Control45_ID_F44_B_cell", "Control46_ID_F47_3_B_cell", "Control47_ID_F55_B_cell", "Case33_ID_F57_B_cell", "Control48_ID_F61_B_cell", "Case34_ID_F62_B_cell", "Case35_ID_F65_B_cell", "Case36_ID_F66_B_cell", "Case37_ID_F67_B_cell", "Case38_ID_F68_B_cell", "Case39_ID_F70_B_cell", "Case40_ID_F72_B_cell", "Case41_ID_F74_B_cell", "Case42_ID_F76_B_cell", "Case43_ID_F77_B_cell", "Case44_ID_F78_B_cell", "Case45_ID_F79_B_cell", "Case46_ID_F80_3_B_cell", "Case47_ID_F81_B_cell", "Case48_ID_F82_B_cell", "Case49_ID_F83_B_cell", "Control49_ID_F86_B_cell", "Control50_ID_F88_B_cell", "Control51_ID_F89_B_cell", "Control52_ID_F90_B_cell", "Control53_ID_F95_B_cell", "Control54_ID_F98_B_cell")
names <- cbind(names, new_names)

pre_name <- names[,2]
new_name <- names[,3]
curr_dir <- getwd()

#Repeat the following code for every flowFrame folder created.

setwd("/Users/orthi/Documents/B_cell_3")
list.files()
file.rename(from = pre_name, to = new_name)
list.files()

setwd(curr_dir)


#PROPORTIONS

x <- gs_pop_get_stats(BMD1_auto_gs)
class(x)

y <- x[946:1386]




#Using ggplot to check if sth wrong with ggcyto

BMD1_matrix <- exprs(BMD1_flowset[[33]])
BMD1_dataframe <- as.data.frame(BMD1_matrix)

ggplot(BMD1_dataframe, aes(x = BMD1_dataframe$`FSC-A`, y = BMD1_dataframe$`SSC-A`)) + geom_hex(bins = 256) + scale_fill_distiller(palette = "Spectral")

#Test negative values in dataset

BMD1_matrix <- exprs(BMD1_flowset[[46]])

each_col(BMD1_flowset[[46]], min)
each_col(BMD1_flowset[[46]], max)

FSCA_vect <- BMD1_matrix[,1]
SSCA_vect <- BMD1_matrix[,4]

length(FSCA_vect[FSCA_vect<0])
length(FSCA_vect[FSCA_vect>0])
length(SSCA_vect[SSCA_vect<0])
length(SSCA_vect[SSCA_vect>0])


#REMOVING NEGATIVE VALUES

install.packages("remotes")
remotes::install_github("ssmpsn2/flowAssist")
library(flowAssist)

BMD1_matrix_46E23 <- exprs(BMD1_flowset[[46]])
BMD1_dataframe_46E23 <- as.data.frame(BMD1_matrix_46E23)

BMD1_df46_pos <- BMD1_dataframe_46E23[BMD1_dataframe_46E23$`FSC-A`>=0, ]
BMD1_df46_allpos <- BMD1_df46_pos[BMD1_df46_pos$`SSC-A`>=0, ]

#DFtoFF is a function in flowAssist for converting dataframe to flowFrame

BMD1_fs46_allpos <- DFtoFF(BMD1_df46_allpos)

ggcyto(BMD1_fs46_allpos, aes(x="FSC-A", y = "SSC-A")) + geom_hex(bins=2000)+ ggcyto_par_set(limits = list(x = c(0, 1e+05), y= c(0, 1e+05)))


#remove doublets

RemoveDoublets(BMD1_compensated, channel1="FSC-A", channel2="FSC-H", nmad=4,
               verbose=FALSE, output="frame")

#If we need to converge multiple flowSets use:

list_1 <- flowSet_to_list(B_cell_gate)
list_2 <- flowSet_to_list(B2_cell_gate)
list_3 <- flowSet_to_list(B3_cell_gate)

large_list <- c(list_1, list_2, list_3)
class(large_list)

BMD_B_cell <- as(large_list, "flowSet")
class(BMD_combined)


#SAVING GATING SETS

save_gs(objectname, path = "where to save")

library(flowWorkspace)
library(ggcyto)
objectname <- load_gs("filepath")


