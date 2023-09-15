#==============================================================================
#    TreeToy_Clone.R : TreeToy Clone
#    Copyright (C) 2023  Bruno Toupance <bruno.toupance@mnhn.fr>
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


# library("RColorBrewer")


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

    return(rtn)
}



#==============================================================================
# check_parameters
#==============================================================================
check_parameters <- function(param_n = 30, param_theta_0 = 10.0, 
    param_growth_factor = 1.0, param_tau = 15.0, max_time = 30)
{

    check_list <- list(check_flag = TRUE, check_msg = "")

    check_list <- check_integer(param_n, label = "n", min_x = 2, max_x = 100, rtn = check_list)
    check_list <- check_real(param_theta_0, label = "Theta_0", min_x = 0, min_inc = TRUE, rtn = check_list)
    check_list <- check_real(param_growth_factor, label = "Growth Factor", min_x = 0, min_inc = TRUE, rtn = check_list)
    check_list <- check_real(param_tau, label = "Tau", min_x = 0, min_inc = TRUE, rtn = check_list)
    check_list <- check_real(max_time, label = "Maximum Time", min_x = 0, min_inc = TRUE, rtn = check_list)

    return(list(msg = check_list$check_msg, flag = check_list$check_flag))
}



#==============================================================================
# INTERNAL FUNCTION - find_position
#==============================================================================
find_position <- function(coal_tree_df, pos_vec, i)
{

    if (!is.null(coal_tree_df)) {

        child_vec <- which(coal_tree_df$parent == i)
        
        if (length(child_vec) == 0) {
            pos_vec[i] <- max(pos_vec) + 1
        } else {
            for (j in 1:length(child_vec)) {
                pos_vec <- find_position(coal_tree_df, pos_vec, child_vec[j])
            }
            pos_vec[i] <- mean(pos_vec[child_vec])
        }
        
        return(pos_vec)
        
    } else {
        return(NULL)
    }
}



#==============================================================================
# branch_color_palette
#==============================================================================
branch_color_palette <- function(coal_tree, branch_color_flag = FALSE)
{
    coal_tree_df <- coal_tree$coal_tree_df
    
    n <- coal_tree$nb_leaf
    nb_col <- n - 1
    
    color_palette_gray <- rep("gray", times = n - 1)
    color_palette_rainbow <- rainbow(n, start = 0, end = 0.7)
    
    color_palette <- color_palette_gray
    if (branch_color_flag) {
        # color_palette <- color_palette_rainbow
        
        # color_palette <- color_palette_gray
        # color_palette_set1 <- brewer.pal(n = 8, name = "Set1")
        # obs <- unique(coal_tree_df$nb_leaf)
        # for (i in 1:min(length(obs), 8)) {
            # color_palette[obs[i]] <- color_palette_set1[i]
        # }

        color_palette <- color_palette_gray
        color_palette[1] <- "red"
    }

    return(color_palette)
}



#==============================================================================
# draw_coal_tree
#==============================================================================
draw_coal_tree <- function(coal_tree, max_time = 3, 
    branch_color_flag = FALSE, add_mutation_flag = TRUE)
{

    if (!is.null(coal_tree)) {
        coal_tree_df <- coal_tree$coal_tree_df

        nb_node <- coal_tree$nb_node
        n <- coal_tree$nb_leaf
        root_i <- coal_tree$root_i

        node_y <- rep(0, times = nb_node)
        node_y <- find_position(coal_tree_df, node_y, root_i)
        node_x <- coal_tree_df$height


        # Create graphic area
        box_y <- c(1, n)
        box_x <- c(0, max_time)
        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))
        plot(box_x, box_y, type = "n", yaxt = "n", bty = "n", main = "", xlab = "", ylab = "")

        branch_col = branch_color_palette(coal_tree, branch_color_flag = branch_color_flag)

        # Plot branches
        for (node_i in 1:(2 * n - 2)) {
            parent_i <- coal_tree_df$parent[node_i]
            segments(node_x[node_i], node_y[node_i], node_x[parent_i], node_y[node_i], 
                col = branch_col[coal_tree_df$nb_leaf[node_i]], lwd = 2)
            segments(node_x[parent_i], node_y[node_i], node_x[parent_i], node_y[parent_i], 
                col = "gray", lwd = 2)
        }

        # Plot root node
        points(node_x[root_i], node_y[root_i], pch = 15, col = "gray", cex = 1)

        # Plot mutations
        if (add_mutation_flag) {
            for (node_i in 1:(2 * n - 2)) {
                parent_i <- coal_tree_df$parent[node_i]
                
                x0 <- node_x[node_i]
                x1 <- node_x[parent_i]
                
                nb_mut <- coal_tree_df$nb_mut[node_i]
                
                if (nb_mut > 0) {
                    
                    x <- coal_tree$mut_height[[node_i]]
                    y <- rep(node_y[node_i], times = nb_mut)
                    
                    # points(x, y, pch = 15, col = "red", cex = 1)
                    points(x, y, pch = 4, col = "black", cex = 1)
                }
            }
        }


        # Add titles
        title(main = "Coalescence Tree")
        title(xlab = "2 x ut")
        text(max_time, n * 0.95, pos = 2, 
            sprintf("TMRCA = %.2f", coal_tree_df$height[root_i]))
    } 
}



#==============================================================================
# compute_pairwise_difference
#==============================================================================
compute_pairwise_difference <- function(coal_tree_df, i, j)
{

    if (!is.null(coal_tree_df)) {

        par_i <- coal_tree_df$parent[i]
        mut_i <- coal_tree_df$nb_mut[i]
        par_j <- coal_tree_df$parent[j]
        mut_j <- coal_tree_df$nb_mut[j]

        while (par_i != par_j) {
            if (coal_tree_df$height[par_i] < coal_tree_df$height[par_j]) {
                mut_i <- mut_i + coal_tree_df$nb_mut[par_i]
                par_i <- coal_tree_df$parent[par_i]
            } else {
                mut_j <- mut_j + coal_tree_df$nb_mut[par_j]
                par_j <- coal_tree_df$parent[par_j]
            }
        }

        return(mut_i + mut_j)
    }
}



#==============================================================================
# compute_mismatch_distribution
#==============================================================================
compute_mismatch_distribution <- function(coal_tree_df)
{
    if (!is.null(coal_tree_df)) {

        nb_node <- nrow(coal_tree_df)
        n <- (nb_node + 1) / 2
        
        nb_diff <- combn(x = n, m = 2, 
            FUN = function(xy){ 
                compute_pairwise_difference(coal_tree_df, xy[1], xy[2])
            }
        )

        return(table(nb_diff))
    }
}



#==============================================================================
# draw_mismatch_distribution
#==============================================================================
draw_mismatch_distribution <- function(coal_tree, max_time = 3, MD_Y_scale_flag = FALSE)
{
    if (!is.null(coal_tree)) {
        coal_tree_df <- coal_tree$coal_tree_df

        nb_node <- coal_tree$nb_node
        n <- coal_tree$nb_leaf

        # Compute the mismatch distribution
        MD_tab <- compute_mismatch_distribution(coal_tree_df)
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

        # Set graphic margins
        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))

        # Plot mismatch distribution
        plot(x, y, type = "h", xlim = c(0, max_time), ylim = c(0, max_y), 
            lwd = 2, col = "gray", 
            main = "", xlab = "", ylab = "")
        points(x, y, pch = 15, col = "gray", cex = 1)

        # Compute the mean of the mismatch distribution
        mean_MD <- sum(x * y)
        abline(v = mean_MD, col = "gray", lty = 2)
        text(max_time, 0.95 * max_y, sprintf("Mean = %.3f", mean_MD), pos = 2)

        # Add titles
        title(main = "Mismatch Distribution")
        title(xlab = "Number of pairwise differences")
        title(ylab = "Frequency")
    }
}



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
# draw_frequency_spectrum
#==============================================================================
draw_frequency_spectrum <- function(coal_tree, 
    DAF_Y_scale_flag = FALSE, DAF_X_scale_flag = FALSE, 
    branch_color_flag = FALSE, std_flag = FALSE)
{

    if (!is.null(coal_tree)) {
        coal_tree_df <- coal_tree$coal_tree_df

        nb_node <- coal_tree$nb_node
        n <- coal_tree$nb_leaf
        
        
        stat_DAF <- rep(coal_tree_df$nb_leaf[-nb_node], coal_tree_df$nb_mut[-nb_node])
        i <- 1:(n - 1)
        ksi_i <- as.numeric(table(factor(stat_DAF, levels = i)))
        stat <- compute_SFS_stat(ksi_i)
        # print(str(stat))

        # Expected DAF
        Ex <- i
        Ey <- 1 / Ex
        Ey <- Ey / sum(Ey)

        if (std_flag) {
            ksi_i <- ksi_i / Ey
            Ey <- Ey / Ey
        }


        # Observed DAF distribution
        x <- i
        y <- ksi_i
        branch_col <- branch_color_palette(coal_tree, branch_color_flag = branch_color_flag)

        # Filter out non-null values
        pos <- which(ksi_i != 0)
        x <- x[pos]
        y <- y[pos]
        branch_col <- branch_col[pos]
        

        # Compute frequencies
        if (stat$stat_S != 0) {
            y <- y / stat$stat_S
        }

        # Set maximum value of DAF
        max_x <- n - 1
        if (DAF_X_scale_flag) {
            max_x <- max(x)
        }

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


        # Set graphic margins
        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))

        # Create graphic area
        if (std_flag) {
            plot(x, y, type = "l", xlim = c(1, max_x), ylim = c(min_y, max_y), 
                lwd = 1, lty = 2, 
                main = "", xlab = "", ylab = "")
        } else {
            plot(x, y, type = "h", xlim = c(1, max_x), ylim = c(min_y, max_y), 
                lwd = 2, col = branch_col, 
                main = "", xlab = "", ylab = "")
        }

        # Add observed DAF distribution
        points(x, y, pch = 15, col = "gray", cex = 1)

        # Add expected DAF distribution
        points(Ex, Ey, pch = 15, col = "black", cex = 1)

        # Add summary statistics values
        text(max_x, 0.95 * max_y, sprintf("S = %d", stat$stat_S), pos = 2)
        # text(max_x, 0.85 * max_y, sprintf("Pi = %.3f", stat$stat_pi), pos = 2)
        text(max_x, 0.85 * max_y, sprintf("Eta = %d", stat$stat_eta), pos = 2)
        text(max_x, 0.75 * max_y, sprintf("Theta_S = %.3f", stat$stat_theta_S), pos = 2)
        text(max_x, 0.65 * max_y, sprintf("Theta_Pi = %.3f", stat$stat_theta_pi), pos = 2)
        text(max_x, 0.55 * max_y, sprintf("Tajima_D = %.3f", stat$stat_tajima_D), pos = 2)

        # Add titles
        title(main = "Derived Allele Frequency Spectrum")
        title(xlab = "DAF (Number of occurences)")
        title(ylab = "Frequency ")

        # No polymorphism message
        if (stat$stat_S == 0) {
            text(0.5 * n, 0.5 * max_y, "No polymorphism", col = "red", adj = 0.5)
        }

    }
}



#==============================================================================
# draw_demography
#==============================================================================
draw_demography <- function(param_theta_0 = 10.0, param_growth_factor = 1.0, 
    param_tau = 15.0, max_time = 30)
{

    # Create graphic area	
    box_x <- c(0, max_time)
    box_y <- c(0, param_theta_0 * max(c(1, param_growth_factor)))
    par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 2, right = 0.5))
    plot(box_x, box_y, type = "n", main = "", xlab = "", ylab = "")

    # Draw demography line 
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

    # Add title
    title(main = "Demography", xlab = "2 x ut", ylab = "Theta")
}



#==============================================================================
# simulate_coal_tree
#==============================================================================
simulate_coal_tree <- function(param_n = 30, param_theta_0 = 10.0, 
    param_growth_factor = 1.0, param_tau = 15.0, 
    ExpTimeFlag = FALSE)
{

    # Parameters cheching
    param_checking <- check_parameters(param_n = param_n, 
        param_theta_0 = param_theta_0, 
        param_growth_factor = param_growth_factor, 
        param_tau = param_tau)
        # print(param_checking)

    if (param_checking$flag) {

        # Number of nodes ("n" external and "n-1" internal)
        nb_node <- 2 * param_n - 1

        if (ExpTimeFlag) {
            coal_time <- param_theta_0 * param_growth_factor / choose(param_n:2, 2)
        } else {
            coal_time <- param_theta_0 * param_growth_factor * rexp(param_n - 1) / choose(param_n:2, 2)
        }

        # Vector storing the height of each node
        height <- c(rep(0, param_n), cumsum(coal_time))

        # Determine the index of nodes occuring before the growth event
        which_pre_growth <- which(height > param_tau)

        # Scale 
        height[which_pre_growth] <- param_tau + (height[which_pre_growth] - param_tau) / param_growth_factor

        # Vector storing the index of the parent of each node
        parent <- rep(NA, nb_node)

        # Vector storing the number of mutations on the branch to the parent of each node 
        nb_mut <- rep(NA, nb_node)

        # Vector storing the number of leaves of each node
        nb_leaf <- c(rep(1, param_n), rep(NA, param_n - 1))

        # Vector storing the indexes of nodes that can coalesce
        avail_node <- 1:param_n 
        for (node_i in (param_n + 1):nb_node) {
            k <- length(avail_node)
            pos <- sample(1:k, size = 2, replace = FALSE)
            parent[avail_node[pos]] <- node_i
            nb_leaf[node_i] <- sum(nb_leaf[avail_node[pos]])
            avail_node <- c(avail_node[-pos], node_i)
        }

        # Add mutations
        mut_height <- vector("list", nb_node)
        for (node_i in 1:(nb_node - 1)) {
            branch_length <- height[parent[node_i]] - height[node_i]
            nb_mut[node_i] <- rpois(1, branch_length / 2)
            mut_height[[node_i]] <- height[node_i] + runif(nb_mut[node_i]) * branch_length
        }


        # Prepare outputs
        coal_tree_df <- data.frame(parent, height, nb_mut, nb_leaf)
        coal_tree_n <- param_n
        coal_tree_nb_node <- nb_node
        coal_tree_root_i <- nb_node
        coal_tree_TMRCA <- height[coal_tree_root_i]
        coal_tree_mut_height <- mut_height

        # Create coalescence tree (R list)
        coal_tree <- list(
          coal_tree_df = coal_tree_df, 
          nb_leaf = coal_tree_n, 
          nb_node = coal_tree_nb_node, 
          root_i = coal_tree_root_i, 
          TMRCA = coal_tree_TMRCA,
          mut_height =  coal_tree_mut_height)

        return(coal_tree)

    } else {
        return(NULL)
    }
}



#==============================================================================
# DoPlot
#==============================================================================
DoPlot <- function(coal_tree, 
    param_n, param_theta_0, param_growth_factor, param_tau, 
    max_time, time_scale_flag = FALSE, MD_Y_scale_flag = FALSE, 
    DAF_X_scale_flag = FALSE, DAF_Y_scale_flag = FALSE, branch_color_flag = FALSE)
{

    # Parameters checking
    param_checking <- check_parameters(
        param_n = param_n, param_theta_0 = param_theta_0, 
        param_growth_factor = param_growth_factor, param_tau = param_tau, 
        max_time = max_time)
        # print(param_checking)

    if (param_checking$flag & !is.null(coal_tree)) {

        coal_tree_df <- coal_tree$coal_tree_df

        # Set graphic layout
        layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))

        # Compute maximum time
        if (time_scale_flag) {
            MD <- compute_mismatch_distribution(coal_tree_df)
            max_MD <- max(as.numeric(names(MD)))
            if (max_MD > 0) {
                max_time <- max_MD
            }
        }

        # Graph #1: Draw coalescence tree
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

        # Graph #2: Draw démography
        draw_demography(param_theta_0 = param_theta_0, param_tau = param_tau, 
            param_growth_factor = param_growth_factor, max_time = max_time)

        # Graph #3: Draw mismatch distribution
        draw_mismatch_distribution(coal_tree, max_time = max_time, 
            MD_Y_scale_flag = MD_Y_scale_flag)

        # Graph #4: Draw allele frequency spectrum
        draw_frequency_spectrum(coal_tree, DAF_Y_scale_flag = DAF_Y_scale_flag, 
            DAF_X_scale_flag = DAF_X_scale_flag, branch_color_flag = branch_color_flag)

        layout(1)

    } else {
        layout(1)

        # Create graphic area
        plot(c(0, 1), c(0, 1), type = "n", 
            xlab = "", ylab = "", main = "", 
            xaxt = "n", yaxt = "n", bty = "n")

        # Build error message
        msg <- sprintf("ERROR: check parameter values%s", param_checking$msg)

        # Display error massage
        text(0.5, 0.5, msg, col = "red", adj = 0.5)

        layout(1)
    }

}



#==============================================================================
# panel_tree_plot
#==============================================================================
panel_tree_plot <- function(coal_tree, 
    param_n, param_theta_0, param_growth_factor, param_tau, 
    max_time, time_scale_flag = FALSE, MD_Y_scale_flag = FALSE, 
    DAF_X_scale_flag = FALSE, DAF_Y_scale_flag = FALSE, branch_color_flag = FALSE)
{

    # Parameters checking
    param_checking <- check_parameters(
        param_n = param_n, param_theta_0 = param_theta_0, 
        param_growth_factor = param_growth_factor, param_tau = param_tau, 
        max_time = max_time)
    # print(param_checking)

    if (param_checking$flag & !is.null(coal_tree)) {

        coal_tree_df <- coal_tree$coal_tree_df

        layout(1)

        max_time <- coal_tree$TMRCA
        draw_coal_tree(coal_tree, max_time = max_time, branch_color_flag = branch_color_flag)

        layout(1)

    } else {
        layout(1)

        # Create graphic area
        plot(c(0, 1), c(0, 1), type = "n", 
            xlab = "", ylab = "", main = "", 
            xaxt = "n", yaxt = "n", bty = "n")

        # Build error message
        msg <- sprintf("ERROR: check parameter values%s", param_checking$msg)

        # Display error massage
        text(0.5, 0.5, msg, col = "red", adj = 0.5)

        layout(1)
    }
}


#==============================================================================
# DoIt
#==============================================================================
DoIt <- function(param_n = 30, param_theta_0 = 10.0, param_growth_factor = 1.0, 
    param_tau = 15.0, max_time = 30, 
    MD_Y_scale_flag = FALSE, DAF_Y_scale_flag = FALSE, time_scale_flag = FALSE, 
    branch_color_flag = FALSE) 
{
    # Generate a coalescence tree with mutations
    coal_tree <- simulate_coal_tree(
        param_n = param_n, param_theta_0 = param_theta_0, 
        param_growth_factor = param_growth_factor, param_tau = param_tau)

    # Plot 4 graphs
    DoPlot(coal_tree, param_n = param_n, param_theta_0 = param_theta_0, 
        param_growth_factor = param_growth_factor, param_tau = param_tau, 
        max_time = max_time, MD_Y_scale_flag = MD_Y_scale_flag, 
        DAF_Y_scale_flag = DAF_Y_scale_flag, 
        time_scale_flag = time_scale_flag, 
        branch_color_flag = branch_color_flag)
}


#==============================================================================

# param_n <- 30; 
# param_theta_0 <- 10.0; 
# param_growth_factor <- 1.0; 
# param_tau <- 15.0; 
# max_time <- 30; 
# MD_Y_scale_flag <- TRUE; 
# DAF_Y_scale_flag <- TRUE; 
# DAF_X_scale_flag <- TRUE; 
# branch_color_flag <- TRUE

# coal_tree <- simulate_coal_tree(param_n = param_n, param_theta_0 = param_theta_0, param_growth_factor = param_growth_factor, param_tau = param_tau)

# layout(1)
# draw_coal_tree(coal_tree, max_time = max_time, branch_color_flag = branch_color_flag)

# layout(1)
# draw_demography(param_theta_0 = param_theta_0, param_tau = param_tau, param_growth_factor = param_growth_factor, max_time = max_time)

# layout(1)
# draw_mismatch_distribution(coal_tree, max_time = max_time, MD_Y_scale_flag = MD_Y_scale_flag)

# layout(1)
# draw_frequency_spectrum(coal_tree, DAF_Y_scale_flag = DAF_Y_scale_flag, DAF_X_scale_flag = DAF_X_scale_flag, branch_color_flag = branch_color_flag)

# DoIt(param_n = param_n, param_theta_0 = param_theta_0, param_growth_factor = param_growth_factor, param_tau = param_tau, max_time = max_time, MD_Y_scale_flag = MD_Y_scale_flag, DAF_Y_scale_flag = DAF_Y_scale_flag, branch_color_flag = branch_color_flag)

