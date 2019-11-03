#==============================================================================
#    ui.R : TreeToy_Clone_Shiny User-Interface
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
				  numericInput(inputId='n',             label='Sample size - integer [min=3 - max=100]:', value=30)
				, numericInput(inputId='Theta0',        label='Theta0 - numeric:', value=10.0)
				, numericInput(inputId='GrowthFactor',  label='GrowthFactor - numeric:', value=1.0)
				, numericInput(inputId='Tau',           label='Tau - numeric:', value=15)
				, actionButton(inputId='go',            label='New Simulation', icon("random"))
			)
			, wellPanel(
				  numericInput(inputId='MaxT',          label='Maximum Time - numeric:', value=50)
				, checkboxInput(inputId='MDScaleFlag',  label='Scale mismatch distribution Y axis', FALSE)
				, checkboxInput(inputId='TimeScaleFlag',label='Scale time X axis', FALSE)
				, checkboxInput(inputId='FSScaleFlag',  label='Scale frequency spectrum Y axis', FALSE)
				, checkboxInput(inputId='DAFScaleFlag', label='Scale DAF X axis', FALSE)
				, checkboxInput(inputId='ColorFlag',    label='Branch colors', FALSE)
			)
		),


#------------------------------------------------------------------------------
# Plot Panel
#------------------------------------------------------------------------------
		mainPanel(
			tabsetPanel(
				  type="tabs"
				, tabPanel(title="MainPlot", plotOutput(outputId="MainPlot", height="600px"))
			)
		)
#------------------------------------------------------------------------------

	)
)

