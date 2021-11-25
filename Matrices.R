#####################################################################################################
##### CODE TO GENERATE FIGURES
#####################################################################################################

### Dependencies

source("Functions.R")

##############################################################################
#### Matrices figures
##############################################################################

ConstC <- BuildC(S = 15, MatType = "Tradeoff", Cd = 2, C0 = 1, bC = 1)
ConstP <- BuildP(S = 15, MatType = "Constant")

CircC <- BuildC(S = 15, MatType = "Circulant", Cd = 2, C0 = 1, bC = 1)
CircP <- BuildP(S = 15, MatType = "Circulant")

RandC <- BuildC(S = 15, MatType = "Random", Cd = 2, C0 = 1, bC = 1)
RandP <- BuildP(S = 15, MatType = "Random")

plConstC <- GridM(ConstC, Title = "(A) Tradeoff\nConsumption", in_color = "Blues")
plConstP <- GridM(ConstP, Title = "(B) Constant\nProduction", in_color = "Greens")
plCircC <- GridM(CircC, Title = "(C) Circulant\nConsumption", in_color = "Blues")
plCircP <- GridM(CircP, Title = "(D) Circulant\nProduction", in_color = "Greens")
plRandC <- GridM(RandC, Title = "C: Random", in_color = "Blues")
plRandP <- GridM(RandP, Title = "P: Random", in_color = "Greens")

grid.arrange(plConstC, plCircC, plRandC, plConstP, plCircP, plRandP, nrow = 2)

tiff("../figs/FigMats.tiff", width = 1200, height = 400, res = 150)
grid.arrange(plConstC, plConstP, plCircC, plCircP, nrow = 1)
dev.off()

### testing consumotion and production structures

BandC <- BuildC(S = 15, MatType = "Banded", Cd = 2, C0 = 1, bC = 1)
ConstP <- BuildP(S = 15, MatType = "Constant")

CorrC <- BuildC(S = 15, MatType = "Correlated", Cd = 2, C0 = 1, bC = 1)
RandP <- BuildP(S = 15, MatType = "Random")

TriC <- BuildC(S = 15, MatType = "Lower Triangular", Cd = 2, C0 = 1, bC = 1)
TriP <- BuildP(S = 15, MatType = "Upper Triangular")

NonSymC <- BuildC(S = 15, MatType = "Non-Sym Tradeoff", Cd = 2, C0 = 1, bC = 1)
SparseC <- BuildC(S = 15, MatType = "Sparse", Cd = 2, C0 = 1, bC = 1)


plBandC <- GridM(BandC, Title = "C: Banded", in_color = "Blues")
plConstP <- GridM(ConstP, Title = "P: Constant", in_color = "Greens")
plCorrC <- GridM(CorrC, Title = "C: Correlated", in_color = "Blues")
plRandP <- GridM(RandP, Title = "P: Random", in_color = "Greens")
plTriC <- GridM(TriC, Title = "C: Lower Triangular", in_color = "Blues")
plTriP <- GridM(TriP, Title = "P: Upper Triangular", in_color = "Greens")
plNonSymC <- GridM(NonSymC, Title = "C: Non-Sym Tradeoff", in_color = "Blues")
plSparseC <- GridM(SparseC, Title = "C: Sparse", in_color = "Blues")


grid.arrange(plBandC, plCorrC, plTriC, plNonSymC, plRandC, plSparseC,
             plConstP, plRandP, plTriP, plConstP, plRandP, plConstP, nrow = 2)


tiff("../figs/SIFigTestMats.tiff", width = 2500, height = 800, res = 150)
grid.arrange(plBandC, plCorrC, plTriC, plNonSymC, plRandC, plSparseC,
             plConstP, plRandP, plTriP, plConstP, plRandP, plConstP, nrow = 2)
dev.off()


