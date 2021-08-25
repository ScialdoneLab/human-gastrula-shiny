library(shinythemes)

ui<-shinyUI(fluidPage(
  theme = shinytheme("flatly"),
  titlePanel(title="Human Gastrulation Data",windowTitle='Human Gastrulation'),
  
  sidebarLayout(
    sidebarPanel(
      
      selectizeInput("var", 
                  label = "Type a gene name and press enter to plot",selected=1,
                 choices=NULL),
      selectizeInput("var2", 
                     label = "Enter another gene to plot in the same UMAP (see Two-Gene expression tab)",selected=1,
                     choices=NULL),
      
      downloadButton("downloadData", label = "Download gene expression matrix(normalized)"),
      downloadButton("downloadUmap", label = "Download annotation and UMAP")
    ),

    
    
    mainPanel(
      navbarPage('',
                  tabPanel("Overview",br(),
                          
                           p("You can access the preprint ", a(href="https://www.biorxiv.org/content/10.1101/2020.07.21.213512v1", "here",target="https://www.biorxiv.org/content/10.1101/2020.07.21.213512v1"),style="text-align:center;color:white",
                           
                         
                             style="text-align:justify;color:white;background-color:rgb(48,62,78);padding:15px;border-radius:15px"),
                           br(),
                           plotOutput("umap_cluster3",width = "85%"),
                           imageOutput("Fig1e"),imageOutput("Fig1a"),imageOutput("Fig1b")
                           ,imageOutput("Fig1c"),imageOutput("Fig1d")),
                           
                  tabPanel("Gene Expression", textOutput("selected_var"),fluid=TRUE,
                           plotOutput("umap_plot",width = "70%"),
                           plotOutput("umap_cluster",width = "85%"),
                           plotOutput("boxplot_cluster")
                          ), 
                  tabPanel("Isoform Expression", 
                           plotlyOutput("iso"),
                           imageOutput("Fig1e_2"),
                           plotOutput("umap_cluster2",width = "85%")
                           
                           ),
                  tabPanel("Two-Gene Expression",  plotOutput("umap_2gene",width = "70%"),
                           plotOutput("umap_2genecolor",width = "40%", height = "200px")),
                  tabPanel("3D diffusion Map",
                           plotlyOutput("diffmap",width="auto")))

      )
      
    )
  )
)


