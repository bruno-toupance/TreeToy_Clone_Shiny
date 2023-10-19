#==============================================================================
#    server.R : TreeToy_Clone_Shiny Server
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


library("shiny")

source("TreeToy_Clone.R")

#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
	function(input, output, session) {

#------------------------------------------------------------------------------
		coal_tree <- reactive(
			{
				tmp <- input$go
				return(
					simulate_coal_tree(
					param_n = input$param_n, 
					param_theta_0 = input$param_theta_0, 
					param_growth_factor = input$param_growth_factor, 
					param_tau = input$param_tau)
				)
			}
		)


#------------------------------------------------------------------------------
		output$MainPlot <- renderPlot(
			{
				Plot <- DoPlot(
					coal_tree(), 
					param_n = input$param_n, 
					param_theta_0 = input$param_theta_0, 
					param_growth_factor = input$param_growth_factor, 
					param_tau = input$param_tau, 
					max_time = input$max_time, 
					time_scale_flag = input$time_scale_flag,
					MD_Y_scale_flag = input$MD_Y_scale_flag, 
					DAF_X_scale_flag = input$DAF_X_scale_flag,
					DAF_Y_scale_flag = input$DAF_Y_scale_flag,
					branch_color_type = input$branch_color_type
					# branch_color_flag = input$branch_color_flag
				)
			}
		)


#------------------------------------------------------------------------------
		output$TreePlot <- renderPlot(
			{
				Plot <- panel_tree_plot(
					coal_tree(), 
					param_n = input$param_n, 
					param_theta_0 = input$param_theta_0, 
					param_growth_factor = input$param_growth_factor, 
					param_tau = input$param_tau, 
					max_time = input$max_time, 
					time_scale_flag = input$time_scale_flag,
					MD_Y_scale_flag = input$MD_Y_scale_flag, 
					DAF_X_scale_flag = input$DAF_X_scale_flag,
					DAF_Y_scale_flag = input$DAF_Y_scale_flag,
					branch_color_type = input$branch_color_type
					# branch_color_flag = input$branch_color_flag
				)
			}
		)


	}
)
#==============================================================================
