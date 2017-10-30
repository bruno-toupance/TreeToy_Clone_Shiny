#==============================================================================
#    ui.R : TreeToy_Clone_Shiny User-Interface
#    Copyright (C) 2017  Bruno Toupance <bruno.toupance@mnhn.fr>
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


#==============================================================================
# shinyUI
#==============================================================================
shinyUI(

	pageWithSidebar(

#------------------------------------------------------------------------------
		headerPanel("TreeToy Clone"),
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Sidebar with input
#------------------------------------------------------------------------------
		sidebarPanel(
			wellPanel(
				  numericInput('n',             'Sample size - integer [3, 100]:', value=30)
				, numericInput('Theta0',        'Theta0 - numeric [0.0, 100.0]:', value=10.0)
				, numericInput('GrowthFactor',  'GrowthFactor - numeric [0.000001, 1000000]:', value=1.0)
				, numericInput('Tau',           'Tau - numeric [0, 1000]:', value=15)
				, numericInput('MaxT',          'Maximum Time - numeric [0, 2000]:', value=50)
				, actionButton('go',            'New Simulation', icon("random"))
			)
			, wellPanel(
				checkboxInput('MDScaleFlag',  'Scale mismatch distribution', FALSE)
				, checkboxInput('FSScaleFlag',  'Scale frequency spectrum', FALSE)
			)
		),


#------------------------------------------------------------------------------
# Plot Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(
				type="tabs"
				, tabPanel("MainPlot", plotOutput("MainPlot", height = "600px"))
			)
		)
#------------------------------------------------------------------------------

	)
)

