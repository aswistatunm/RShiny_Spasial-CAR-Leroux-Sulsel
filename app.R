#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(rgdal)
library(spdep)
library(CARBayes)
library(readxl)
library(leaflet)
library(DT) #Tampilan data Tabel
library(coda) #uji konvergensi menggunakan trace plot

# ui Objek
#Tampilan Halaman
#dashboardHeader=Judul Kolom
#dashboardSidebar=Judul Baris

# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader (title =  " SPATIAL ANALYSIS "),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName="data"),
      menuItem(" Model Comparison", tabName = "model"),
      menuItem("Spatial Model", tabName = "model1"),
      menuItem("Convergence Test", tabName = "plot"),
      menuItem("Map of Relative Risk", tabName="peta")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName="data",
              fluidRow(
                box( fileInput(
                  inputId = "filedata1",
                  label = "Upload Data. Reading Data From Excel File",
                  accept = c(".xlsx")))
              ),
              fluidRow(
                DTOutput("table"), title = "Data")
      ),
      tabItem(tabName="model",
              fluidRow(
                box(width = 6,
                    selectInput("data1", "Cases",
                                choices=c(), selected =c())),
                box(width = 6,
                    selectInput("data2", "Population",
                                choices=c(), selected =c()))
              ),
              fluidRow(titlePanel("Model Comparison"), DTOutput("modelx"))
      ),
      tabItem(tabName = "model1",
              fluidRow(titlePanel("Hyperprior"),
                       box(width = 3,numericInput("prior",NULL,0.1, min = 0, max = 1)),
                       box(width = 3,numericInput("tau",NULL,0.1, min = 0, max = 1))
              ),
              fluidRow(DTOutput("mo"),title = "Parameter Estimation"),
              fluidRow(DTOutput("mo1"), title = "Model fit")
      ),
      tabItem(tabName = "plot",
              fluidRow( titlePanel("Trace Plot for tau2 Samples"),
                        plotOutput("plot1")),
              fluidRow( titlePanel("Trac Plot for rho Samples"),
                        plotOutput("plot2")),
              fluidRow( titlePanel("Trace Plot for beta Samples"),
                        plotOutput("plot3")),
              fluidRow( titlePanel("Trace Plot for phi Samples"),
                        plotOutput("plot4"))
      ),
      
      tabItem(tabName="peta",
              fluidRow( box(width = 12, height = 430,leafletOutput("petak"))
              ),
              fluidRow(
                DTOutput("table1"), title = "Data Resiko Relatif")
      )
    ))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # ===== Input =====
  
  Data <- reactive({
    req(input$filedata1)
    read_excel (input$filedata1$datapath)
  })
  observeEvent(Data(), {
    updateSelectInput(session, "data1", 
                      choices = names(Data()), selected = names(Data()))
  })
  observeEvent(Data(), {
    updateSelectInput(session, "data2", 
                      choices = names(Data()), selected = names(Data()))
  })
  
  Map <- readOGR("SulSel/sul-sel.shp")
  
  names(Map@data)[1] <-"Districts"
  
  # ===== Proses =====
  
  W.nb <- poly2nb(Map)
  W <- nb2mat(W.nb, style = "B")
  W
  
  
  
  Covid <- reactive({
    Covid <- Data()
    Exp<-sum(Covid[input$data1])*Covid[input$data2]/sum(Covid[input$data2])
    Covid$Exp <- Exp[1:dim(Covid)[1],]
    Covid <- as.data.frame(Covid)
    Covid
  })
  
  # Define user-defined function Residuals MMI
  get.morans.I.modified <- function(x, h){
    
    Wm <- t(h / colSums(h))
    
    # Compute modified moran's I
    Z<-(x-mean(x))
    ZW <- Wm %*% x-mean(x)
    a <- t(Z) %*% ZW
    b <- (t(Z) %*% Z)^(1/2)
    c <- (t(ZW) %*% ZW)^(1/2)
    
    I.modified <- (a/(b*c))
    return(I.modified)
  }
  
  c = matrix(c(1,1,0.1,0.5,0.5,0.01,0.1,0.1,0.05,0.0005), nrow = 5, byrow = FALSE)
  v = matrix(0,nrow = 5, ncol = 3)
  colnames(v) = c("DIC","WAIC","MMI")
  row.names(v) = c("Leroux Model IG(1, 0.01)","Leroux Model IG(1, 0.1)",
                   "Leroux Model IG(0.1, 0.1)","Leroux Model IG(0.5, 0.05)",
                   "Leroux Model IG(0.5, 0.0005)")
  
  G <- reactive({
    Formula <-Covid()[input$data1][1:dim(Covid())[1],]~offset(log(Covid()$Exp))
    for (i in 1:dim(c)[1]) {
      f = as.vector(c[i,])
      set.seed(1)
      ModelLeroux <-S.CARleroux(formula=Formula,
                                data=Covid(),
                                family="poisson",
                                W=W,
                                burnin=10000,
                                prior.tau2=f,
                                n.sample=110000)
      
      x <-ModelLeroux$residuals[,1]
      
      v[i,1] <- as.matrix(round(ModelLeroux$modelfit,2))[1,1]
      v[i,2] <- as.matrix(round(ModelLeroux$modelfit,2))[3,1]
      v[i,3] <- as.matrix(round(get.morans.I.modified(x, W),2))
    }
    v
  })
  
  
  ModelLeroux <- reactive({
    Formula <-Covid()[input$data1][1:dim(Covid())[1],]~offset(log(Covid()$Exp))
    set.seed(1)
    ModelLeroux <-S.CARleroux(formula=Formula,
                              data=Covid(),
                              family="poisson",
                              W=W,
                              burnin=10000,
                              prior.tau2=c(input$prior,input$tau),
                              n.sample=110000)
    
    ModelLeroux   
  })
  
  par(mfrow=c(1,2),mar=c(4,4,1.5,1.5),mex=0.8,tcl=0.3)
  
  tau2.samples<- reactive({
    mcmc.list(ModelLeroux()$samples$tau2)
  })
  
  rho.samples<- reactive({
    mcmc.list(ModelLeroux()$samples$rho)
  })
  
  beta.samples<- reactive({
    mcmc.list(ModelLeroux()$samples$beta)
  })
  phi.samples<- reactive({
    mcmc.list(ModelLeroux()$samples$phi)
  })  
  
  # Relative Risk
  MapS <- reactive({
    Covid1 <- Covid()
    Covid1$RR<-round(ModelLeroux()$fitted.values/Covid1$Exp,2)
    Map <- merge(Map,Covid1)
    Map@data <- Map@data[c(1,3,12,13,15)]
    Map
  }) 
  
  
  
  
  
  
  
  # ===== Output =====
  output$table <- renderDT(Data())
  
  output$modelx <- renderDT(as.data.frame(G()))
  
  output$mo <- renderDT(as.data.frame(ModelLeroux()$summary.results))
  
  output$mo1 <- renderDT({
    Modelfit <- round(ModelLeroux()$modelfit,2)
    t(as.data.frame( Modelfit))
  })
  
  output$plot1 <- renderPlot( plot(tau2.samples()))
  
  output$plot2 <- renderPlot( plot(rho.samples()))
  
  output$plot3 <- renderPlot( plot(beta.samples()) )
  
  output$plot4 <- renderPlot( plot(phi.samples()) )
  
  output$table1 <- renderDT(MapS()@data)
  
  output$petak <- renderLeaflet({
    if (is.null(Data()) | is.null(MapS())) {
      return(NULL)
    }
    Map1 <- MapS()
    
    # Create variableplot
    Map1$variableplot <- as.numeric(Map1@data$RR)
    
    # Create leaflet
    pal <- colorBin(colorRamp(c("Blue","White","Red")),domain = Map1$variableplot, bins = 6)
    labels <- sprintf("%s: %g",Map1@data$Districts , Map1$variableplot) %>%
      lapply(htmltools::HTML)
    
    leaflet(Map1) %>%
      addTiles() %>%
      addPolygons(
        fillColor = ~ pal(variableplot),
        color = "white",
        dashArray = "3",
        fillOpacity = 0.7,
        label = labels
      ) %>%
      leaflet::addLegend(
        pal = pal, values = ~variableplot,
        opacity = 0.7, title = NULL
      )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
