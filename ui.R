#==============================================================================
#    ui.R: TreeToy_Clone_Shiny User-Interface
#    Copyright (C) 2025  Bruno Toupance <bruno.toupance@mnhn.fr>
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
        
        # Header ---
        headerPanel("TreeToy Clone"),
        
        
        # Sidebar with input ---
        sidebarPanel(
            wellPanel(
                numericInput(
                    inputId = 'param_n', 
                    label = 'Sample size - integer [min = 2 - max = 100]:', 
                    value = 30),
                # bsTooltip(
                #     id = "param_n",
                #     title = "Integer [min = 2 - max = 100]"),
                
                numericInput(
                    inputId = 'param_theta_0', 
                    label = 'Theta0 - numeric:', value = 20.0),
                
                numericInput(
                    inputId = 'param_growth_factor', 
                    label = 'Growth Factor - numeric:', value = 1.0),
                
                numericInput(
                    inputId = 'param_tau', 
                    label = 'Tau - numeric:', value = 15),
                
                actionButton(
                    inputId = 'go', 
                    label = 'New Simulation', icon("random"))
            ),
            
            wellPanel(
                numericInput(
                    inputId = 'max_time', 
                    label = 'Maximum Time - numeric:', value = 50),
                
                checkboxInput(
                    inputId = 'time_scale_flag',
                    label = 'Scale time X axis', value = FALSE),
                
                checkboxInput(
                    inputId = 'MD_Y_scale_flag', 
                    label = 'Scale mismatch distribution Y axis', 
                    value = TRUE),
                
                checkboxInput(
                    inputId = 'DAF_X_scale_flag', 
                    label = 'Scale DAF X axis', value = FALSE),
                
                checkboxInput(
                    inputId = 'DAF_Y_scale_flag', 
                    label = 'Scale DAF Y axis', value = TRUE),
                
                checkboxInput(
                    inputId = 'Expected_WF_DAF_flag', 
                    label = 'Show expected WF DAF', value = FALSE),
                
                radioButtons(
                    inputId = "branch_color_type", 
                    label = "Branch colors:", 
                    choices = c("none", "all", "singleton"))
                
            )
        ),
        
        
        # Plot Panel ---
        mainPanel(
            tabsetPanel(
                type = "tabs",
                
                tabPanel(
                    title = "Main", 
                    plotOutput(outputId = "main_plot", height = "600px")
                ),
                
                tabPanel(
                    title = "Tree", 
                    plotOutput(outputId = "tree_plot", height = "600px")
                ),
                
            )
        )
        
    )
)

