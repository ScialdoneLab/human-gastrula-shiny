options(shiny.sanitize.errors = F) 
source("./global.R")

server<-function(input,output,session){
  updateSelectizeInput(session=session, 'var', choices =  c(Choose = '', gene.names), server = TRUE,selected="TBXT")
  updateSelectizeInput(session=session, 'var2', choices =  c(Choose = '', gene.names), server = TRUE,selected="MESP1")
  
  output$umap_plot <- renderPlot({
    validate(
      need(input$var, "Please select another gene")
    )
    df_scatter<-data.table(x=umap$X0,
                           y=umap$X1,
                           z=unlist(express_vals[[input$var]]))
    plot.std.col(df_scatter, xname="", yname="", title=paste0(input$var))
  })
  
  output$iso <- renderPlotly({
    validate(
      need(input$var, "Please select another gene")
    )
    #is_df$cluster <-is_df$cluster
    
    plot_iso(input$var)
    
  })
  
  
  output$umap_cluster3<-output$umap_cluster2<-output$umap_cluster <- renderPlot({
    plot_cluster()
    
  })
  output$boxplot_cluster<-renderPlot({
    df_boxplot<-data.table(z=unlist(express_vals[[input$var]]),
                           cluster=umap$cluster_id)
    boxplot_cluster(df_boxplot,input$var)
    
    
  })
  output$umap_2gene<-renderPlot({
    plot_2gene(input$var,input$var2)
  })
  
  output$umap_2genecolor<-renderPlot({
    plot_2genecolor(input$var,input$var2)
  })
  
    
  output$diffmap<-renderPlotly({
    plot_3d()
    
    
  })
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("expression_values", "rds", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste0(base,'express_vals.rds'), file)
    }
  )
  output$downloadUmap <- downloadHandler(
    filename <- function() {
      paste("annot_umap", "rds", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste0(base,'umap.rds'), file)
    }
  )
  
  
  output$Fig1e_2<-output$Fig1e <- renderImage({
    filename <- 'www/spat.png'
    list(src = filename,
         width = 600,height = 350)
    
  }, deleteFile = FALSE)
  output$Fig1a <- renderImage({
    filename <- 'www/Figure_1_DissectionRegions.png'
    list(src = filename,
         width = 600,height = 350)
  }, deleteFile = FALSE)
  output$Fig1b <- renderImage({
    filename <- 'www/Figure_1_WholeEmbryo.png'
    list(src = filename,
         width = 600,height = 350)
   
  }, deleteFile = FALSE)
  output$Fig1c <- renderImage({
    filename <- 'www/Figure_1_DorsalView.png'
    list(src = filename,
         width = 600,height = 350)
    
  }, deleteFile = FALSE)
  output$Fig1d <- renderImage({
    filename <- 'www/Figure_1_VentralView.png'
    list(src = filename,
         width = 600,height = 350,align='center')
    
  }, deleteFile = FALSE)
  
}


