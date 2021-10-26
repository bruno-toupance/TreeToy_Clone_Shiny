#==============================================================================
#    server.R : TreeToy_Clone_Shiny Server
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


library(shiny)

source("TreeToy_Clone.R")

#==============================================================================
# shinyServer
#==============================================================================
shinyServer(
	function(input, output, session) {

#------------------------------------------------------------------------------
		CoalTree <- reactive({
			Tmp <- input$go
			return(simulate_coal_tree(
				n = input$n, 
				Theta0 = input$Theta0, 
				GrowthFactor = input$GrowthFactor, 
				Tau = input$Tau))
		})
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
		output$MainPlot <- renderPlot({
			Plot <- DoPlot(CoalTree(), 
				n = input$n, 
				Theta0 = input$Theta0, 
				GrowthFactor = input$GrowthFactor, 
				Tau = input$Tau, 
				MaxT = input$MaxT, 
				MDScaleFlag = input$MDScaleFlag, 
				TimeScaleFlag = input$TimeScaleFlag,
				FSScaleFlag = input$FSScaleFlag,
				DAFScaleFlag = input$DAFScaleFlag,
				ColorFlag = input$ColorFlag)
		})
#------------------------------------------------------------------------------

	}
)
#==============================================================================
