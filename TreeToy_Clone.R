#==============================================================================
#    TreeToy_Clone.R : TreeToy Clone
#    Copyright (C) 2019  Bruno Toupance <bruno.toupance@mnhn.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================



#==============================================================================
# INTERNAL FUNCTION - find_position
#==============================================================================
find_position <- function(CoalTree, Pos, i)
{
#------------------------------------------------------------------------------
	Child <- which(CoalTree$Parent == i)
	if (length(Child) == 0) {
		Pos[i] <- max(Pos)+1
	} else {
		for (j in 1:length(Child)) {
			Pos <- find_position(CoalTree, Pos, Child[j])
		}
		Pos[i] <- mean(Pos[Child])
	}
	return(Pos)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_coal_tree
#==============================================================================
draw_coal_tree <- function(CoalTree, MaxT=3, ColorFlag=FALSE)
{
#------------------------------------------------------------------------------
	m <- nrow(CoalTree)
	n <- (m+1)/2
#------------------------------------------------------------------------------
	y <- rep(0, m)
	y <- find_position(CoalTree, y, m)
	x <- CoalTree$Height
#------------------------------------------------------------------------------
	if (ColorFlag) {
		Color <- rainbow(n-1, start=0, end=0.7)
	} else {
		Color <- rep("gray", n-1)
	}
#------------------------------------------------------------------------------
	BoxY <- c(1, n)
	BoxX <- c(0, MaxT)
	par(mar=c(bottom=0.5+4,left=0.5+4,top=0.5+2,right=0.5))
	plot(BoxX, BoxY, type="n", yaxt="n", bty="n", main="", xlab="", ylab="")
#------------------------------------------------------------------------------
	for (i in 1:(2*n-2)) {
		Parent <- CoalTree$Parent[i]
		segments(x[i], y[i], x[Parent], y[i], col=Color[CoalTree$NbLeaf[i]], 
			lwd=2)
		segments(x[Parent], y[i], x[Parent], y[Parent], col="gray", lwd=2)
	}
#------------------------------------------------------------------------------
	points(x[m], y[m], pch=15, col="gray", cex=1)  # TMRCA
#------------------------------------------------------------------------------
	title(main="Coalescence Tree")
	title(xlab="2ut")
	text(MaxT, n*0.95, sprintf("TMRCA = %.2f", CoalTree$Height[m]), pos=2)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# compute_pairwise_difference
#==============================================================================
compute_pairwise_difference <- function(CoalTree, i, j)
{
#------------------------------------------------------------------------------
	p1 <- CoalTree$Parent[i]; m1 <- CoalTree$Mutation[i]
	p2 <- CoalTree$Parent[j]; m2 <- CoalTree$Mutation[j]
	while (p1 != p2) {
		if (CoalTree$Height[p1] < CoalTree$Height[p2]) {
			m1 <- m1+CoalTree$Mutation[p1]; p1 <- CoalTree$Parent[p1]
		} else {
			m2 <- m2+CoalTree$Mutation[p2]; p2 <- CoalTree$Parent[p2]
		}
	}
#------------------------------------------------------------------------------
	return(m1+m2)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# compute_mismatch_distribution
#==============================================================================
compute_mismatch_distribution <- function(CoalTree)
{
#------------------------------------------------------------------------------
	m <- nrow(CoalTree)
	n <- (m+1)/2
	NbDiff <- combn(n, m=2, FUN=function(pair){ 
		compute_pairwise_difference(CoalTree, pair[1], pair[2]) } )
#------------------------------------------------------------------------------
	return(table(NbDiff))
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_mismatch_distribution
#==============================================================================
draw_mismatch_distribution <- function(CoalTree, MaxT=3, MDScaleFlag=FALSE)
{
#------------------------------------------------------------------------------
	m <- nrow(CoalTree)
	n <- (m+1)/2
	Tab <- compute_mismatch_distribution(CoalTree)
	Tab <- prop.table(Tab)
	x <- as.numeric(names(Tab))
	y <- as.vector(Tab)
	MaxY <- 1
	if (MDScaleFlag) {
		MaxY <- max(y)
		if (MaxY == 0) {
			MaxY <- 1
		}
	}
	par(mar=c(bottom=0.5+4,left=0.5+4,top=0.5+2,right=0.5))
	plot(x, y, type="h", xlim=c(0, MaxT), ylim=c(0, MaxY), lwd=2, col="gray", 
		main="", xlab="", ylab="")
	points(x, y, pch=15, col="gray", cex=1)
	MeanMD <- sum(x*y)
	abline(v=MeanMD, col="gray", lty=2)
	text(MaxT, 0.95*MaxY, sprintf("Mean = %.3f", MeanMD), pos=2)
	title(main="Mismatch Distribution")
	title(xlab="Number of pairwise differences")
	title(ylab="Frequency")
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# compute_Tajima_D
#==============================================================================
compute_Tajima_D <- function(n, Pi, S)
{
#------------------------------------------------------------------------------
	i <- 1:(n-1)
	a1 <- sum(1/i)
	a2 <- sum(1/i^2)
	b1 <- (n+1)/(3*(n-1))
	b2 <- (2*(n^2+n+3))/(9*n*(n-1))
	c1 <- b1-1/a1
	c2 <- b2-(n+2)/(a1*n)+a2/a1^2
	e1 <- c1/a1
	e2 <- (c2)/(a1^2+a2)
	TajD <- (Pi-S/a1)/sqrt(e1*S+e2*S*(S-1))
#------------------------------------------------------------------------------
	return(TajD)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# compute_SFS_stat
#==============================================================================
compute_SFS_stat <- function(Ksi_i)
{
#------------------------------------------------------------------------------
	n <- length(Ksi_i)+1
	i <- 1:(n-1)
	a1 <- sum(1/i)
	a2 <- sum(1/i^2)
	S <- sum(Ksi_i)
	ThetaS <- S/a1
	Pi <- sum(Ksi_i*i*(n-i))/choose(n, 2)
	ThetaPi <- Pi
	Eta <- Ksi_i[1]
	TajD <- compute_Tajima_D(n, Pi, S)
#------------------------------------------------------------------------------
	return(list(S=S, Pi=Pi, Eta=Eta, 
		ThetaS=ThetaS, ThetaPi=ThetaPi, 
		TajD=TajD))
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_frequency_spectrum
#==============================================================================
draw_frequency_spectrum <- function(CoalTree, FSScaleFlag=FALSE, 
	DAFScaleFlag=FALSE, ColorFlag=FALSE)
{
#------------------------------------------------------------------------------
	m <- nrow(CoalTree)
	n <- (m+1)/2
	DAF <- rep(CoalTree$NbLeaf[-m], CoalTree$Mutation[-m])
	i <- 1:(n-1)
	Ksi_i <- as.numeric(table(factor(DAF, levels=i)))
	Stat <- compute_SFS_stat(Ksi_i)
	# print(str(Stat))
#------------------------------------------------------------------------------
	Ex <- i
	Ey <- 1/Ex; Ey <- Ey/sum(Ey)
#------------------------------------------------------------------------------
	if (ColorFlag) {
		Color <- rainbow(n-1, start=0, end=0.7)
	} else {
		Color <- rep("gray", n-1)
	}
#------------------------------------------------------------------------------
	x <- i
	y <- Ksi_i
	x <- x[(Ksi_i != 0)]
	y <- y[(Ksi_i != 0)]
#------------------------------------------------------------------------------
	Color <- Color[(Ksi_i != 0)]
#------------------------------------------------------------------------------
	if (Stat$S != 0) {
		y <- y/Stat$S
	}
#------------------------------------------------------------------------------
	MaxX <- n-1
	if (DAFScaleFlag) {
		MaxX <- max(x)
	}
#------------------------------------------------------------------------------
	MaxY <- 1
	if (FSScaleFlag) {
		MaxY <- max(c(y, Ey))
	}
#------------------------------------------------------------------------------
	par(mar=c(bottom=0.5+4,left=0.5+4,top=0.5+2,right=0.5))
	plot(x, y, type="h", xlim=c(1, MaxX), ylim=c(0, MaxY), lwd=2, col=Color, 
		main="", xlab="", ylab="")
	points(x, y, pch=15, col="gray", cex=1)
#------------------------------------------------------------------------------
	points(Ex, Ey, pch=15, col="black", cex=1)  # Expected SFS
#------------------------------------------------------------------------------
	text(MaxX, 0.95*MaxY, sprintf("S = %d", Stat$S), pos=2)
	# text(MaxX, 0.85*MaxY, sprintf("Pi = %.3f", Stat$Pi), pos=2)
	text(MaxX, 0.85*MaxY, sprintf("Eta = %d", Stat$Eta), pos=2)
	text(MaxX, 0.75*MaxY, sprintf("Theta_S = %.3f", Stat$ThetaS), pos=2)
	text(MaxX, 0.65*MaxY, sprintf("Theta_Pi = %.3f", Stat$ThetaPi), pos=2)
	text(MaxX, 0.55*MaxY, sprintf("Tajima_D = %.3f", Stat$TajD), pos=2)
	title(main="Derived Allele Frequency Spectrum")
	title(xlab="DAF (Number of occurences)")
	title(ylab="Frequency ")
#------------------------------------------------------------------------------
	if (Stat$S == 0) {
		text(0.5*n, 0.5*MaxY, "No polymorphism", col="red", adj=0.5)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_demography
#==============================================================================
draw_demography <- function(Theta0=10.0, GrowthFactor=1.0, Tau=15.0, MaxT=30)
{
#------------------------------------------------------------------------------
	BoxX <- c(0, MaxT)
	BoxY <- c(0, Theta0*max(c(1, GrowthFactor)))
	par(mar=c(bottom=0.5+4,left=0.5+4,top=0.5+2,right=0.5))
	plot(BoxX, BoxY, type="n", main="", xlab="", ylab="")
#------------------------------------------------------------------------------
	if (GrowthFactor > 1.0) {
		lines(c(0, Tau, Tau, MaxT), Theta0*c(GrowthFactor, GrowthFactor, 1, 1), 
			lty=1, col="red", lwd=3)
	} else {
		if (GrowthFactor < 1.0) {
			lines(c(0, Tau, Tau, MaxT), Theta0*c(GrowthFactor, 
				GrowthFactor, 1, 1), lty=1, col="blue", lwd=3)
		} else {
			lines(c(0, MaxT), c(Theta0, Theta0), lty=1, col="black", lwd=3)
		}
	}
#------------------------------------------------------------------------------
	title(main="Demography", xlab="2ut", ylab="Theta")
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# check_parameters
#==============================================================================
check_parameters <- function(n=30, Theta0=10.0, GrowthFactor=1.0, Tau=15.0, 
	MaxT=30)
{
#------------------------------------------------------------------------------
	Flag <- TRUE
	Msg <- ""
#------------------------------------------------------------------------------
	if (is.numeric(n)) {
		if (n == round(n)) {
			if ( (n < 3) | (n > 100) ) {
				Flag <- FALSE
				Msg <- sprintf("%s\n%s", Msg, "FAIL: [n] out of bounds")
			}
		} else {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [n] not integer")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [n] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(Theta0)) {
		if ( (Theta0 <= 0) ) {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [Theta0] out of bounds")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [Theta0] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(GrowthFactor)) {
		if ( (GrowthFactor <= 0) ) {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [GrowthFactor] out of bounds")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [GrowthFactor] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(Tau)) {
		if ( (Tau <= 0) ) {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [Tau] out of bounds")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [Tau] not numeric")
	}
#------------------------------------------------------------------------------
	if (is.numeric(MaxT)) {
		if ( (MaxT <= 0) ) {
			Flag <- FALSE
			Msg <- sprintf("%s\n%s", Msg, "FAIL: [MaxT] out of bounds")
		}
	} else {
		Flag <- FALSE
		Msg <- sprintf("%s\n%s", Msg, "FAIL: [MaxT] not numeric")
	}
#------------------------------------------------------------------------------
	return(list(Msg=Msg, Flag=Flag))
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# simulate_coal_tree
#==============================================================================
simulate_coal_tree <- function(n=30, Theta0=10.0, GrowthFactor=1.0, Tau=15.0, 
	ExpTimeFlag=FALSE)
{
#------------------------------------------------------------------------------
	ParamChecking <- check_parameters(n=n, Theta0=Theta0, 
		GrowthFactor=GrowthFactor, Tau=Tau)
	# print(ParamChecking)
#------------------------------------------------------------------------------
	if (ParamChecking$Flag) {
		m <- 2*n-1
		if (ExpTimeFlag) {
			CoalTime <- Theta0*GrowthFactor/choose(n:2, 2)
		} else {
			CoalTime <- Theta0*GrowthFactor*rexp(n-1)/choose(n:2, 2)
		}
		Height <- c(rep(0, n), cumsum(CoalTime))
		PreGrowth <- which(Height > Tau)
		Height[PreGrowth] <- Tau + (Height[PreGrowth]-Tau)/GrowthFactor
		Parent <- rep(NA, m)
		Mutation <- rep(NA, m)
		NbLeaf <- c(rep(1, n), rep(NA, n-1))
		Available <- 1:n 
		for (i in (n+1):m) {
			k <- length(Available)
			Pos <- sample(1:k, size=2, replace=FALSE)
			Parent[Available[Pos]] <- i
			NbLeaf[i] <- sum(NbLeaf[Available[Pos]])
			Available <- c(Available[-Pos], i)
		}
		for (i in 1:(m-1)) {
			BranchLength <- Height[Parent[i]] - Height[i]
			Mutation[i] <- rpois(1, BranchLength/2)
		}
		return(data.frame(Parent, Height, Mutation, NbLeaf))
#------------------------------------------------------------------------------
	} else {
		return(NULL)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# DoPlot
#==============================================================================
DoPlot <- function(CoalTree, n, Theta0, GrowthFactor, Tau, MaxT, 
	MDScaleFlag=FALSE, TimeScaleFlag=FALSE, FSScaleFlag=FALSE, 
	DAFScaleFlag=FALSE, ColorFlag=FALSE)
{
#------------------------------------------------------------------------------
	ParamChecking <- check_parameters(n=n, Theta0=Theta0, 
		GrowthFactor=GrowthFactor, Tau=Tau, MaxT=MaxT)
	# print(ParamChecking)
#------------------------------------------------------------------------------
	if (ParamChecking$Flag & !is.null(CoalTree)) {
		layout(matrix(1:4, nrow=2, ncol=2, byrow=TRUE))
		
		if (TimeScaleFlag) {
			MD <- compute_mismatch_distribution(CoalTree)
			MaxMD <- max(as.numeric(names(MD)))
			if (MaxMD > 0) {
				MaxT = MaxMD
			}
		}
		draw_coal_tree(CoalTree, MaxT=MaxT, ColorFlag=ColorFlag)
		if (GrowthFactor > 1.0) {
			abline(v=Tau, col="red", lty=2)
		} else {
			if (GrowthFactor < 1.0) {
				abline(v=Tau, col="blue", lty=2)
			} else {
				E_TMRCA <- 2*Theta0*(1-1/n)
				abline(v=E_TMRCA, col="black", lty=2)
				text(MaxT, n*0.85, sprintf("E[TMRCA] = %.2f", E_TMRCA), pos=2)
			}
		}
		
		draw_demography(Theta0=Theta0, Tau=Tau, GrowthFactor=GrowthFactor, 
			MaxT=MaxT)
		draw_mismatch_distribution(CoalTree, MaxT=MaxT, 
			MDScaleFlag=MDScaleFlag)
		draw_frequency_spectrum(CoalTree, FSScaleFlag=FSScaleFlag, 
			DAFScaleFlag=DAFScaleFlag, ColorFlag=ColorFlag)
	
		layout(1)
#------------------------------------------------------------------------------
	} else {
		plot(c(0, 1), c(0, 1), type="n", xlab="", ylab="", main="", 
			xaxt="n", yaxt="n", bty="n")
		Msg <- sprintf("ERROR: check parameter values%s", ParamChecking$Msg)
		text(0.5, 0.5, Msg, col="red", adj=0.5)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# DoIt
#==============================================================================
DoIt <- function(n=30, Theta0=10.0, GrowthFactor=1.0, Tau=15.0, MaxT=30, 
	MDScaleFlag=FALSE, FSScaleFlag=FALSE, TimeScaleFlag=FALSE, 
	ColorFlag=FALSE) 
{
#------------------------------------------------------------------------------
	CoalTree <- simulate_coal_tree(n=n, Theta0=Theta0, 
		GrowthFactor=GrowthFactor, Tau=Tau)
	DoPlot(CoalTree, n=n, Theta0=Theta0, GrowthFactor=GrowthFactor, Tau=Tau, 
		MaxT=MaxT, MDScaleFlag=MDScaleFlag, FSScaleFlag=FSScaleFlag, 
		TimeScaleFlag=TimeScaleFlag, ColorFlag=ColorFlag)
#------------------------------------------------------------------------------
}
#==============================================================================


# n <- 30; Theta0 <- 10.0; GrowthFactor <- 1.0; Tau <- 15.0; MaxT <- 30; MDScaleFlag <- TRUE; FSScaleFlag <- TRUE; ColorFlag <- TRUE
# DoIt(n=n, Theta0=Theta0, GrowthFactor=GrowthFactor, Tau=Tau, MaxT=MaxT, MDScaleFlag=MDScaleFlag, FSScaleFlag=FSScaleFlag, ColorFlag=ColorFlag)
