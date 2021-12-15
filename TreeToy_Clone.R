#==============================================================================
#    TreeToy_Clone.R : TreeToy Clone
#    Copyright (C) 2021  Bruno Toupance <bruno.toupance@mnhn.fr>
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
find_position <- function(coal_tree, pos_vec, i)
{
#------------------------------------------------------------------------------
	child_vec <- which(coal_tree$parent == i)
	if (length(child_vec) == 0) {
		pos_vec[i] <- max(pos_vec) + 1
	} else {
		for (j in 1:length(child_vec)) {
			pos_vec <- find_position(coal_tree, pos_vec, child_vec[j])
		}
		pos_vec[i] <- mean(pos_vec[child_vec])
	}
	return(pos_vec)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_coal_tree
#==============================================================================
draw_coal_tree <- function(coal_tree, max_time = 3, branch_color_flag = FALSE)
{
#------------------------------------------------------------------------------
	nb_node <- nrow(coal_tree)
	n <- (nb_node + 1) / 2
	root_i <- nb_node
#------------------------------------------------------------------------------
	y <- rep(0, times = nb_node)
	y <- find_position(coal_tree, y, root_i)
	x <- coal_tree$height
#------------------------------------------------------------------------------
	if (branch_color_flag) {
		branch_col <- rainbow(n - 1, start = 0, end = 0.7)
	} else {
		branch_col <- rep("gray", n - 1)
	}
#------------------------------------------------------------------------------
	box_y <- c(1, n)
	box_x <- c(0, max_time)
	par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))
	plot(box_x, box_y, type = "n", yaxt = "n", bty = "n", main = "", xlab = "", ylab = "")
#------------------------------------------------------------------------------
	for (i in 1:(2 * n - 2)) {
		parent <- coal_tree$parent[i]
		segments(x[i], y[i], x[parent], y[i], 
			col = branch_col[coal_tree$nb_leaf[i]], lwd = 2)
		segments(x[parent], y[i], x[parent], y[parent], 
			col = "gray", lwd = 2)
	}
#------------------------------------------------------------------------------
	points(x[root_i], y[root_i], pch = 15, col = "gray", cex = 1)  # TMRCA
#------------------------------------------------------------------------------
	title(main = "Coalescence Tree")
	title(xlab = "2 x ut")
	text(max_time, n * 0.95, pos = 2, 
		sprintf("TMRCA = %.2f", coal_tree$height[root_i]))
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# compute_pairwise_difference
#==============================================================================
compute_pairwise_difference <- function(coal_tree, i, j)
{
#------------------------------------------------------------------------------
	par_i <- coal_tree$parent[i]
	mut_i <- coal_tree$nb_mut[i]
	par_j <- coal_tree$parent[j]
	mut_j <- coal_tree$nb_mut[j]
	while (par_i != par_j) {
		if (coal_tree$height[par_i] < coal_tree$height[par_j]) {
			mut_i <- mut_i + coal_tree$nb_mut[par_i]
			par_i <- coal_tree$parent[par_i]
		} else {
			mut_j <- mut_j + coal_tree$nb_mut[par_j]
			par_j <- coal_tree$parent[par_j]
		}
	}
#------------------------------------------------------------------------------
	return(mut_i + mut_j)
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# compute_mismatch_distribution
#==============================================================================
compute_mismatch_distribution <- function(coal_tree)
{
#------------------------------------------------------------------------------
	nb_node <- nrow(coal_tree)
	n <- (nb_node + 1) / 2
	nb_diff <- combn(x = n, m = 2, 
			FUN = function(xy){ 
				compute_pairwise_difference(coal_tree, xy[1], xy[2])
			}
		)
#------------------------------------------------------------------------------
	return(table(nb_diff))
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_mismatch_distribution
#==============================================================================
draw_mismatch_distribution <- function(coal_tree, max_time = 3, MD_Y_scale_flag = FALSE)
{
	nb_node <- nrow(coal_tree)
	n <- (nb_node + 1) / 2
	MD_tab <- compute_mismatch_distribution(coal_tree)
	MD_tab <- prop.table(MD_tab)
	x <- as.numeric(names(MD_tab))
	y <- as.vector(MD_tab)
	max_y <- 1
	if (MD_Y_scale_flag) {
		max_y <- max(y)
		if (max_y == 0) {
			max_y <- 1
		}
	}
	par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))

	plot(x, y, type = "h", xlim = c(0, max_time), ylim = c(0, max_y), 
		lwd = 2, col = "gray", 
		main = "", xlab = "", ylab = "")
	points(x, y, pch = 15, col = "gray", cex = 1)
	
	mean_MD <- sum(x * y)
	abline(v = mean_MD, col = "gray", lty = 2)
	text(max_time, 0.95 * max_y, sprintf("Mean = %.3f", mean_MD), pos = 2)
	
	title(main = "Mismatch Distribution")
	title(xlab = "Number of pairwise differences")
	title(ylab = "Frequency")
}
#==============================================================================



#==============================================================================
# compute_Tajima_D
#==============================================================================
compute_Tajima_D <- function(n, stat_pi, stat_S)
{
	i <- 1:(n - 1)
	a1 <- sum(1 / i)
	a2 <- sum(1 / i^2)
	b1 <- (n + 1) / (3 * (n - 1))
	b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
	c1 <- b1 - 1 / a1
	c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
	e1 <- c1 / a1
	e2 <- (c2) / (a1^2 + a2)
	stat_tajima_D <- (stat_pi - stat_S / a1) / sqrt(e1 * stat_S + e2 * stat_S * (stat_S - 1))

	return(stat_tajima_D)
}
#==============================================================================



#==============================================================================
# compute_SFS_stat
#==============================================================================
compute_SFS_stat <- function(ksi_i)
{
	n <- length(ksi_i) + 1
	i <- 1:(n - 1)
	a1 <- sum(1 / i)
	a2 <- sum(1 / i^2)
	stat_S <- sum(ksi_i)
	stat_theta_S <- stat_S / a1
	stat_pi <- sum(ksi_i * i * (n - i)) / choose(n, 2)
	stat_theta_pi <- stat_pi
	stat_eta <- ksi_i[1]
	stat_tajima_D <- compute_Tajima_D(n, stat_pi, stat_S)

	return(list(stat_S = stat_S, stat_pi = stat_pi, stat_eta = stat_eta, 
		stat_theta_S = stat_theta_S, stat_theta_pi = stat_theta_pi, 
		stat_tajima_D = stat_tajima_D))
}
#==============================================================================



#==============================================================================
# draw_frequency_spectrum
#==============================================================================
draw_frequency_spectrum <- function(coal_tree, 
	DAF_Y_scale_flag = FALSE, DAF_X_scale_flag = FALSE, 
	branch_color_flag = FALSE, std_flag = FALSE)
{
#------------------------------------------------------------------------------
	nb_node <- nrow(coal_tree)
	n <- (nb_node + 1) / 2
	stat_DAF <- rep(coal_tree$nb_leaf[-nb_node], coal_tree$nb_mut[-nb_node])
	i <- 1:(n - 1)
	ksi_i <- as.numeric(table(factor(stat_DAF, levels = i)))
	stat <- compute_SFS_stat(ksi_i)
	# print(str(stat))
#------------------------------------------------------------------------------
	Ex <- i
	Ey <- 1 / Ex
	Ey <- Ey / sum(Ey)
#------------------------------------------------------------------------------
	if (std_flag) {
		ksi_i <- ksi_i / Ey
		Ey <- Ey / Ey
	}
#------------------------------------------------------------------------------
	if (branch_color_flag) {
		branch_col <- rainbow(n - 1, start = 0, end = 0.7)
	} else {
		branch_col <- rep("gray", n - 1)
	}
#------------------------------------------------------------------------------
	x <- i
	y <- ksi_i
	x <- x[(ksi_i != 0)]
	y <- y[(ksi_i != 0)]
#------------------------------------------------------------------------------
	branch_col <- branch_col[(ksi_i != 0)]
#------------------------------------------------------------------------------
	if (stat$stat_S != 0) {
		y <- y / stat$stat_S
	}
#------------------------------------------------------------------------------
	max_x <- n - 1
	if (DAF_X_scale_flag) {
		max_x <- max(x)
	}
#------------------------------------------------------------------------------
	if (std_flag) {
		min_y <- min(y)
		max_y <- max(y)
	} else {
		min_y <- 0
		max_y <- 1
		if (DAF_Y_scale_flag) {
			max_y <- max(c(y, Ey))
		}
	}
#------------------------------------------------------------------------------
	par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))
	if (std_flag) {
		plot(x, y, type = "l", xlim = c(1, max_x), ylim = c(min_y, max_y), 
			lwd = 1, lty = 2, 
			main = "", xlab = "", ylab = "")
	} else {
		plot(x, y, type = "h", xlim = c(1, max_x), ylim = c(min_y, max_y), 
			lwd = 2, col = branch_col, 
			main = "", xlab = "", ylab = "")
	}
	points(x, y, pch = 15, col = "gray", cex = 1)
#------------------------------------------------------------------------------
	points(Ex, Ey, pch = 15, col = "black", cex = 1)  # Expected SFS
#------------------------------------------------------------------------------
	text(max_x, 0.95 * max_y, sprintf("S = %d", stat$stat_S), pos = 2)
	# text(max_x, 0.85 * max_y, sprintf("Pi = %.3f", stat$stat_pi), pos = 2)
	text(max_x, 0.85 * max_y, sprintf("Eta = %d", stat$stat_eta), pos = 2)
	text(max_x, 0.75 * max_y, sprintf("Theta_S = %.3f", stat$stat_theta_S), pos = 2)
	text(max_x, 0.65 * max_y, sprintf("Theta_Pi = %.3f", stat$stat_theta_pi), pos = 2)
	text(max_x, 0.55 * max_y, sprintf("Tajima_D = %.3f", stat$stat_tajima_D), pos = 2)
	title(main = "Derived Allele Frequency Spectrum")
	title(xlab = "DAF (Number of occurences)")
	title(ylab = "Frequency ")
#------------------------------------------------------------------------------
	if (stat$stat_S == 0) {
		text(0.5 * n, 0.5 * max_y, "No polymorphism", col = "red", adj = 0.5)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# draw_demography
#==============================================================================
draw_demography <- function(param_theta_0 = 10.0, param_growth_factor = 1.0, 
	param_tau = 15.0, max_time = 30)
{
#------------------------------------------------------------------------------
	box_x <- c(0, max_time)
	box_y <- c(0, param_theta_0 * max(c(1, param_growth_factor)))
	par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))
	plot(box_x, box_y, type = "n", main = "", xlab = "", ylab = "")
#------------------------------------------------------------------------------
	if (param_growth_factor > 1.0) {
		lines(c(0, param_tau, param_tau, max_time), 
			param_theta_0 * c(param_growth_factor, param_growth_factor, 1, 1), 
			lty = 1, col = "red", lwd = 3)
	} else {
		if (param_growth_factor < 1.0) {
			lines(c(0, param_tau, param_tau, max_time), 
				param_theta_0 * c(param_growth_factor, param_growth_factor, 1, 1),
				lty = 1, col = "blue", lwd = 3)
		} else {
			lines(c(0, max_time), c(param_theta_0, param_theta_0), lty = 1, col = "black", lwd = 3)
		}
	}
#------------------------------------------------------------------------------
	title(main = "Demography", xlab = "2 x ut", ylab = "Theta")
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# check_integer
#==============================================================================
check_integer <- function(
	x = 0,
	label = "x", 
	min_x = -Inf, 
	max_x = +Inf, 
	min_inc = FALSE,
	max_inc = FALSE,
	rtn = list(check_flag = TRUE, check_msg = "")) 
{
	if (is.numeric(x)) {
		if (x == round(x)) {
			out_of_bounds_flag <- FALSE
			if (min_inc) {
				if (x <= min_x) {
					out_of_bounds_flag <- TRUE
				}
			} else {
				if (x < min_x) {
					out_of_bounds_flag <- TRUE
				}
			}
			if (max_inc) {
				if (x >= max_x) {
					out_of_bounds_flag <- TRUE
				}
			} else {
				if (x > max_x) {
					out_of_bounds_flag <- TRUE
				}
			}
			if (out_of_bounds_flag) {
				rtn$check_flag <- FALSE
				rtn$check_msg <- sprintf("%s\nFAIL: [%s] out of bounds", rtn$check_msg, label)
			}
		} else {
			rtn$check_flag <- FALSE
			rtn$check_msg <- sprintf("%s\nFAIL: [%s] not integer", rtn$check_msg, label)
		}
	} else {
		rtn$check_flag <- FALSE
		rtn$check_msg <- sprintf("%s\nFAIL: [%s] not numeric", rtn$check_msg, label)
	}
#------------------------------------------------------------------------------
	return(rtn)
}


#==============================================================================
# check_real
#==============================================================================
check_real <- function(
	x = 0.0,
	label = "x", 
	min_x = -Inf, 
	max_x = +Inf, 
	min_inc = FALSE,
	max_inc = FALSE,
	rtn = list(check_flag = TRUE, check_msg = "")) 
{
	if (is.numeric(x)) {
		out_of_bounds_flag <- FALSE
		if (min_inc) {
			if (x <= min_x) {
				out_of_bounds_flag <- TRUE
			}
		} else {
			if (x < min_x) {
				out_of_bounds_flag <- TRUE
			}
		}
		if (max_inc) {
			if (x >= max_x) {
				out_of_bounds_flag <- TRUE
			}
		} else {
			if (x > max_x) {
				out_of_bounds_flag <- TRUE
			}
		}
		if (out_of_bounds_flag) {
			rtn$check_flag <- FALSE
			rtn$check_msg <- sprintf("%s\nFAIL: [%s] out of bounds", rtn$check_msg, label)
		}
	} else {
		rtn$check_flag <- FALSE
		rtn$check_msg <- sprintf("%s\nFAIL: [%s] not numeric", rtn$check_msg, label)
	}
#------------------------------------------------------------------------------
	return(rtn)
}



#==============================================================================
# check_parameters
#==============================================================================
check_parameters <- function(param_n = 30, param_theta_0 = 10.0, 
	param_growth_factor = 1.0, param_tau = 15.0, max_time = 30)
{
#------------------------------------------------------------------------------
	check_list <- list(check_flag = TRUE, check_msg = "")
#------------------------------------------------------------------------------
	check_list <- check_integer(param_n, label = "n", min_x = 2, max_x = 100, rtn = check_list)
	check_list <- check_real(param_theta_0, label = "Theta_0", min_x = 0, min_inc = TRUE, rtn = check_list)
	check_list <- check_real(param_growth_factor, label = "Growth Factor", min_x = 0, min_inc = TRUE, rtn = check_list)
	check_list <- check_real(param_tau, label = "Tau", min_x = 0, min_inc = TRUE, rtn = check_list)
	check_list <- check_real(max_time, label = "Maximum Time", min_x = 0, min_inc = TRUE, rtn = check_list)
#------------------------------------------------------------------------------
	return(list(msg = check_list$check_msg, flag = check_list$check_flag))
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# simulate_coal_tree
#==============================================================================
simulate_coal_tree <- function(param_n = 30, param_theta_0 = 10.0, 
	param_growth_factor = 1.0, param_tau = 15.0, 
	ExpTimeFlag = FALSE)
{
#------------------------------------------------------------------------------
	param_checking <- check_parameters(param_n = param_n, param_theta_0 = param_theta_0, 
		param_growth_factor = param_growth_factor, param_tau = param_tau)
	# print(param_checking)
#------------------------------------------------------------------------------
	if (param_checking$flag) {
		nb_node <- 2 * param_n - 1
		if (ExpTimeFlag) {
			CoalTime <- param_theta_0 * param_growth_factor / choose(param_n:2, 2)
		} else {
			CoalTime <- param_theta_0 * param_growth_factor * rexp(param_n - 1) / choose(param_n:2, 2)
		}
		height <- c(rep(0, param_n), cumsum(CoalTime))
		PreGrowth <- which(height > param_tau)
		height[PreGrowth] <- param_tau + (height[PreGrowth] - param_tau) / param_growth_factor
		parent <- rep(NA, nb_node)
		nb_mut <- rep(NA, nb_node)
		nb_leaf <- c(rep(1, param_n), rep(NA, param_n - 1))
		avail_node <- 1:param_n 
		for (i in (param_n + 1):nb_node) {
			k <- length(avail_node)
			pos <- sample(1:k, size = 2, replace = FALSE)
			parent[avail_node[pos]] <- i
			nb_leaf[i] <- sum(nb_leaf[avail_node[pos]])
			avail_node <- c(avail_node[-pos], i)
		}
		for (i in 1:(nb_node - 1)) {
			BranchLength <- height[parent[i]] - height[i]
			nb_mut[i] <- rpois(1, BranchLength / 2)
		}
		return(data.frame(parent, height, nb_mut, nb_leaf))
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
DoPlot <- function(coal_tree, 
	param_n, param_theta_0, param_growth_factor, param_tau, 
	max_time, time_scale_flag = FALSE, MD_Y_scale_flag = FALSE, 
	DAF_X_scale_flag = FALSE, DAF_Y_scale_flag = FALSE, branch_color_flag = FALSE)
{
#------------------------------------------------------------------------------
	param_checking <- check_parameters(
		param_n = param_n, param_theta_0 = param_theta_0, 
		param_growth_factor = param_growth_factor, param_tau = param_tau, 
		max_time = max_time)
	# print(param_checking)
#------------------------------------------------------------------------------
	if (param_checking$flag & !is.null(coal_tree)) {
		layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
		
		if (time_scale_flag) {
			MD <- compute_mismatch_distribution(coal_tree)
			max_MD <- max(as.numeric(names(MD)))
			if (max_MD > 0) {
				max_time = max_MD
			}
		}
		draw_coal_tree(coal_tree, max_time = max_time, branch_color_flag = branch_color_flag)
		if (param_growth_factor > 1.0) {
			abline(v = param_tau, col = "red", lty = 2)
		} else {
			if (param_growth_factor < 1.0) {
				abline(v = param_tau, col = "blue", lty = 2)
			} else {
				E_TMRCA <- 2 * param_theta_0 * (1 - 1 / param_n)
				abline(v = E_TMRCA, col = "black", lty = 2)
				text(max_time, param_n * 0.85, sprintf("E[TMRCA] = %.2f", E_TMRCA), pos = 2)
			}
		}
		
		draw_demography(param_theta_0 = param_theta_0, param_tau = param_tau, 
			param_growth_factor = param_growth_factor, max_time = max_time)
			
		draw_mismatch_distribution(coal_tree, max_time = max_time, 
			MD_Y_scale_flag = MD_Y_scale_flag)
			
		draw_frequency_spectrum(coal_tree, DAF_Y_scale_flag = DAF_Y_scale_flag, 
			DAF_X_scale_flag = DAF_X_scale_flag, branch_color_flag = branch_color_flag)
	
		layout(1)
#------------------------------------------------------------------------------
	} else {
		plot(c(0, 1), c(0, 1), type = "n", 
			xlab = "", ylab = "", main = "", 
			xaxt = "n", yaxt = "n", bty = "n")
		msg <- sprintf("ERROR: check parameter values%s", param_checking$msg)
		text(0.5, 0.5, msg, col = "red", adj = 0.5)
	}
#------------------------------------------------------------------------------
}
#==============================================================================



#==============================================================================
# DoIt
#==============================================================================
DoIt <- function(param_n = 30, param_theta_0 = 10.0, param_growth_factor = 1.0, 
	param_tau = 15.0, max_time = 30, 
	MD_Y_scale_flag = FALSE, DAF_Y_scale_flag = FALSE, time_scale_flag = FALSE, 
	branch_color_flag = FALSE) 
{
#------------------------------------------------------------------------------
	coal_tree <- simulate_coal_tree(
		param_n = param_n, param_theta_0 = param_theta_0, 
		param_growth_factor = param_growth_factor, param_tau = param_tau)
		
	DoPlot(coal_tree, param_n = param_n, param_theta_0 = param_theta_0, 
		param_growth_factor = param_growth_factor, param_tau = param_tau, 
		max_time = max_time, MD_Y_scale_flag = MD_Y_scale_flag, 
		DAF_Y_scale_flag = DAF_Y_scale_flag, 
		time_scale_flag = time_scale_flag, 
		branch_color_flag = branch_color_flag)
#------------------------------------------------------------------------------
}
#==============================================================================


# param_n <- 30; param_theta_0 <- 10.0; param_growth_factor <- 1.0; param_tau <- 15.0; max_time <- 30; MD_Y_scale_flag <- TRUE; DAF_Y_scale_flag <- TRUE; branch_color_flag <- TRUE
# coal_tree <- simulate_coal_tree(param_n = param_n, param_theta_0 = param_theta_0, param_growth_factor = param_growth_factor, param_tau = param_tau)
# draw_coal_tree(coal_tree, max_time = max_time, branch_color_flag = branch_color_flag)
# draw_demography(param_theta_0 = param_theta_0, param_tau = param_tau, param_growth_factor = param_growth_factor, max_time = max_time)
# draw_mismatch_distribution(coal_tree, max_time = max_time, MD_Y_scale_flag = MD_Y_scale_flag)
# draw_frequency_spectrum(coal_tree, DAF_Y_scale_flag = DAF_Y_scale_flag, DAF_X_scale_flag = DAF_X_scale_flag, branch_color_flag = branch_color_flag)
# DoIt(param_n = param_n, param_theta_0 = param_theta_0, param_growth_factor = param_growth_factor, param_tau = param_tau, max_time = max_time, MD_Y_scale_flag = MD_Y_scale_flag, DAF_Y_scale_flag = DAF_Y_scale_flag, branch_color_flag = branch_color_flag)
