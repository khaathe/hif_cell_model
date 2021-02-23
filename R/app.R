## app.R ##
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(R.matlab)
library(ggplot2)
#library(sp)
library(tidyr)
library(dplyr)
#library(plot3D)
library(viridis)
#library(oce)
library(reshape2)
library(xml2)
library(plotly)

### DATA INIT
sliceZ <- 0
#setwd("../")
#getwd()
print(getwd())
folder <- "../output/"
##GET FIRST XML
data_param <- read_xml(paste0(folder, "output00000000.xml"))
measured_data <-
  data_param %>% xml_find_all("//label") %>% xml_text()
envir_data <-
  data_param %>% xml_find_all("//variable") %>% xml_attr(attr = "name")

## RETRIEVE MEDIUM'S AND CELL'S PARAMETERS NAMES
cell_measured_data_name = c()
for (i in 1:length(measured_data)) {
  j = 1
  vector_name = c(".x", ".y", ".z")
  vector_size = (data_param %>% xml_find_all("//label") %>% xml_attr(attr =
                                                                       "size") %>% as.numeric())[i]
  while (j <= vector_size) {
    if (vector_size > 1) {
      cell_measured_data_name <-
        c(cell_measured_data_name,
          paste0(measured_data[i], vector_name[j]))
    }
    else{
      cell_measured_data_name <-
        c(cell_measured_data_name, measured_data[i])
    }
    j = j + 1
  }
}



## R PLOT FUNCTIONS

GET_TIMESTEP <- function() {
  return (as.numeric(gsub(".xml", "", gsub(
    "output", "", dir(folder, pattern = "*.xml")
  ))))
}



GET_MEDIUM_OBJECT <- function(TIME) {
  files.names <- dir(folder, pattern = "*_microenvironment0.mat")
  path_file <-
    paste(folder, files.names[TIME], sep = "")
  MEDIUM <- readMat(path_file)
  MEDIUM <- MEDIUM$multiscale.microenvironment
  MEDIUM <- t(MEDIUM)
  MEDIUM <- data.frame(MEDIUM, TIME)
  colnames(MEDIUM) <-
    c("x", "y", "z", "constant", envir_data, "TIME")
  if ("H" %in% colnames(MEDIUM)) {
    MEDIUM$pH <- -log10(MEDIUM$H / 1000)
  }
  MEDIUM$Distance <-
    sqrt((MEDIUM$x - 0) ^ 2 + (MEDIUM$y - 0) ^ 2 + (MEDIUM$z - 0) ^
           2 * 1.0)
  
  return(MEDIUM)
}


GET_CELL_OBJECT <- function(TIME) {
  files.names2 <- dir(folder, pattern = "*_cells_physicell.mat")
  cells <- data.frame()
  path_file <-
    paste(folder, files.names2[TIME], sep = "")
  cells_temp <- readMat(path_file)
  cells_temp <- cells_temp$cells
  cells_temp <- t(cells_temp)
  cells_temp <- as.data.frame(cells_temp)
  colnames(cells_temp) <- cell_measured_data_name
  cells_temp$distance <-
    with(cells_temp, sqrt((cells_temp$position.x ^ 2) + (position.y ^ 2) + (position.z ^ 2)))
  cells_temp$simulation_time <- TIME
  cells <- rbind(cells, cells_temp)
  return(cells)
}



EnvironmentMAP <- function(PARAMETER, TIME) {
  MEDIUM_OBJECT <- GET_MEDIUM_OBJECT(TIME)
  rangePARAMETER <- range(MEDIUM_OBJECT[, PARAMETER])
  scalePARAMETER <-
    scale_fill_viridis_c(PARAMETER,
                         option = "plasma",
                         limits = c(floor(rangePARAMETER[1]), ceiling(rangePARAMETER[2])))
  themeNoLegends <-
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1, 0)
    )
  
  ggplot(data = MEDIUM_OBJECT, aes(x = x, y = y, fill = MEDIUM_OBJECT[, PARAMETER])) + geom_raster(interpolate = TRUE) +
    scalePARAMETER + themeNoLegends + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    coord_equal()
  # path <-
  #   paste("./output/pngtemp/",PARAMETER,"/", TIME, ".png", sep = "")
  # print(path)
  # ggsave(
  #   path,
  #   width = 10,
  #   height = 10,
  #   dpi = 72,
  #   units = "in",
  #   device = "png"
  # )
}



plotCellatTime <- function(PARAMETER, time) {
  CELL_OBJECT <- GET_CELL_OBJECT(time)
  CELL_OBJECT <-
    CELL_OBJECT[, c(
      "position.x",
      "position.y",
      "position.z",
      "distance",
      PARAMETER,
      "simulation_time"
    )]
  plot <- ggplot(data = CELL_OBJECT,
                 aes_string(x = "distance", y = PARAMETER)) +
    geom_point(col = "red", alpha = 0.5) +
    geom_smooth(alpha = 0.02) +
    ggtitle(PARAMETER, "/Distance") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 13)) +
    xlab("Distance") + ylab(PARAMETER) + theme(axis.text.x = element_text(angle = 0, hjust = 1))
  return(plot)
}

plotSpheroidAtTIME <- function(PARAMETER, time) {
  CELL_OBJECT <- GET_CELL_OBJECT(time)
  
  size_of_cells <-
    0.01 * (0.5 * (3 * CELL_OBJECT$total_volume / (4 * pi)) ^ (1 / 3))
  
  #return(
    # plot_ly(
    #   CELL_OBJECT,
    #   x = ~ position.x,
    #   y = ~ position.y,
    #   size = size_of_cells,
    #   type="scatter",
    #   mode="markers",
    #   color = CELL_OBJECT[[PARAMETER]],
    #   marker = list(
    #     symbol = 'circle',
    #     opacity = 1,
    #     sizemode = "diameter"
    #   ),
    #   sizes = c(0.5, 4))%>%
    # layout(autosize=F, width=500, height=500))
  
  ggplot(CELL_OBJECT,
        aes(position.x, position.y)) + geom_point(aes(colour = CELL_OBJECT[[PARAMETER]]),
                                                  alpha = 0.8,
                                                  size = 1.5 * (0.5 * (3 * CELL_OBJECT$total_volume / (4 * pi)) ^
                                                                  (1 / 3))) + scale_colour_viridis(option = "plasma",
                                                                                                  direction = -1,
                                                                                                  discrete = FALSE) + coord_equal() + labs(colour = "")
}

plotSpheroid3DAtTIME <- function(PARAMETER, time) {
  CELL_OBJECT <- GET_CELL_OBJECT(time)
  size_of_cells <-
    0.01 * (0.5 * (3 * CELL_OBJECT$total_volume / (4 * pi)) ^ (1 / 3))
  
  return(
    plot_ly(
      CELL_OBJECT,
      x = ~ position.x,
      y = ~ position.y,
      z = ~ position.z,
      size = size_of_cells,
      color = CELL_OBJECT[[PARAMETER]],
      marker = list(
        symbol = 'circle',
        opacity = 1,
        sizemode = "diameter"
      ),
      sizes = c(0.5, 4)
    )
  )
  
}






















### UI
ui <- dashboardPage(
  dashboardHeader(title = "PhysiCell Simulations"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      style = "position:fixed;width: 220px;",
      menuItem(
        "Dashboard",
        tabName = "dashboard",
        icon = icon("dashboard")
      ),
      #menuItem("Widgets", tabName = "widgets", icon = icon("th"))
      #sliderInput("slider_time", "Timestep:", min=min(GET_TIMESTEP()), max=max(GET_TIMESTEP()), step=1, value=max(GET_TIMESTEP()))
      div("TIMELINE", style = "text-align:center"),
      checkboxInput("keep_on_last_step", "Keep at last step", FALSE),
      noUiSliderInput(
        inputId = "slider_time",
        label = NULL,
        min = min(GET_TIMESTEP()),
        max = max(GET_TIMESTEP()),
        step = 1,
        value = max(GET_TIMESTEP()),
        margin = 0,
        orientation = "vertical",
        direction = "rtl",
        width = "300px",
        height = "200px"
      )
    )
  ),
  ## Body content
  dashboardBody(tabItems(
    # First tab content
    tabItem(tabName = "dashboard",
            fluidPage(
              column(6,
                     fluidRow(
                       box(
                         selectInput("Input_p1", "Parameters", envir_data, multiple = FALSE),
                         plotOutput("plot1"),
                         width = 6
                       ),
                       box(
                         selectInput("Input_p2", "Parameters", measured_data, multiple = FALSE),
                         plotOutput("plot2"),
                         width = 6
                       )),
                       fluidRow(
                       box(
                         selectInput("Input_p3", "Parameters", measured_data, multiple = FALSE),
                         plotOutput("plot3"),
                         width = 6
                       ),
                       box(
                         selectInput("Input_p5", "Parameters", measured_data, multiple = FALSE),
                         plotOutput("plot5"),
                         width = 6
                       )
                     )),
                     column(6, 
                            fluidRow(
                            box(
                         selectInput("Input_p4", "Parameters", measured_data, multiple = FALSE),
                         plotlyOutput("plot4", height = "900px"),
                         width = 12,
                         height="100%"
                       )
                     )
                     ),
              column(12,
                     fluidRow(box(
                       selectInput("Input_p6", "Parameters", measured_data, multiple = FALSE),
                       plotOutput("plot6"),
                       width = 12
                     ))))
            ),
    # Second tab content
    tabItem(tabName = "widgets",
            h2("Widgets tab content"))
  ))
)





server <- function(input, output, session) {
  autoInvalidate <- reactiveTimer(1000, session)
  observeEvent({
    autoInvalidate()
    input$keep_on_last_step
    1
  }, {
    if (input$keep_on_last_step == FALSE) {
      updateNoUiSliderInput(session,
                            "slider_time",
                            value = NULL,
                            range = c(0, max(GET_TIMESTEP())))
    }
    else{
      updateNoUiSliderInput(
        session,
        "slider_time",
        value = max(GET_TIMESTEP()),
        range = c(0, max(GET_TIMESTEP()))
      )
    }
    
  })
  
  
  output$plot1 <- renderPlot({
    EnvironmentMAP("oxygen", input$slider_time + 1)
  })
  output$plot2 <- renderPlot({
    plotSpheroidAtTIME(input$Input_p2, input$slider_time + 1)
  })
  output$plot3 <- renderPlot({
    plotCellatTime(input$Input_p3, input$slider_time + 1)
  })
  
  output$plot4 <- renderPlotly({
    plotSpheroid3DAtTIME(input$Input_p4, input$slider_time + 1)
  })
  
  output$plot5 <- renderPlot({
    plotSpheroidAtTIME(input$Input_p5, input$slider_time + 1)
  })
  
  output$plot6 <- renderPlot({
    plotSpheroidAtTIME(input$Input_p6, input$slider_time + 1)
  })
  
  
}



shinyApp(ui, server)