pplotmod<-function(df=df){
    library(ggplot2)
    library(plotly)
    # Ejemplo:
    #df <- data.frame(
    #  designation= paste0("gen_",1:100),
    #  repF=rep(1:2,50),
    #   Observed = rnorm(100, 10, 2),
    #   fitted = rnorm(100, 10, 2),
    #   residuals = rnorm(100),
    #  outlier=rep(c("GOOD","outlier"),50)
    # )
      #df$residuals=df$stdt
    #=======================================
    # 1. residuals vs fitted
    #=======================================
    p <- ggplot(df, aes(x = fitted,y = residuals)) +
      geom_point(aes(color= outlier,
          text = paste("ID:", idRow,"<br>designation:", designation,"<br>rep:", repF,"<br>fitted:", round(fitted, 2),"<br>residuals:", round(residuals, 2))
        ),alpha = 0.7) + 
      scale_color_manual(
        values = c( "GOOD" = "#A4D65E","outlier" = "#FF7F7F"))+
      geom_smooth(method = "loess",se = FALSE,color = "black")+#0072B2") +
      geom_hline(yintercept = 0,linetype = "dashed",color = "red") +
      theme_bw()
      #labs(title = "residuals vs Fitted",x = "Fitted values",y = "residuals")
    g1<-ggplotly(p, tooltip = "text")
    #=======================================
    # 2. QQ Plot
    #=======================================
    qq <- qqnorm(df$stdt, plot.it = FALSE)
    qq_df <- data.frame(theoretical = qq$x,sample = qq$y,designation = df$designation,residuals = df$stdt,outlier=df$outlier)
    
    p2 <- ggplot(qq_df,aes(x = theoretical,y = sample)) +
      geom_point(aes(color=outlier,text = paste("designation:", designation,"<br>residuals:", round(residuals, 3)))) +
      scale_color_manual(values = c( "GOOD" = "#A4D65E","outlier" = "#FF7F7F"))+
      geom_abline(slope = 1,intercept = 0,color = "red") +
      theme_bw() 
      #labs(title = "QQ Plot of Residuals",x = "Theoretical Quantiles",y = "Sample Quantiles")
    g2<-ggplotly(p2, tooltip = "text")%>% style(showlegend = FALSE)
    #=======================================
    # 3. Histogram of Residuals
    #=======================================
    p3 <- ggplot(df,aes(x = stdt,text = paste("residuals:", round(stdt,2)))) +
      geom_histogram(bins = 50) +
      theme_bw()
    #  labs(title = "Histogram of Residuals",x = "residuals",y = "Frequency") +
      
    g3<-ggplotly(p3, tooltip = "text")
    #=======================================
    # 4. Observed vs fitted
    #=======================================
    #p5 <- ggplot(df,aes(x = fitted,y = Observed,text = paste("Observed:", round(Observed,2),"<br>Fitted:", round(fitted,2)))) +
    #  geom_point(size = 2, aes(color=outlier))+
    #  scale_color_manual(values = c( "GOOD" = "#A4D65E","outlier" = "#FF7F7F"))+
    #  geom_abline(slope = 1,intercept = 0,linetype = "dashed",color = "red") +
    #  theme_bw()
    #g5<-ggplotly(p5, tooltip = "text")%>% style(showlegend = FALSE)
    
    #=======================================
    # 5. Residuals vs Leverage
    #=======================================
    pn<-unique(df$p)
    N<-unique(df$N)
    lev_lim<-(2*pn)/N
    cook_lim=4/N
    lev_max=max(max(df$leverage),lev_lim)
    h <- seq(0.001,(lev_max+0.2),length.out = N)
    cook_curve <- function(D, p, h){
      sqrt(D * p * (1 - h) / h)
    }
    c05 <- cook_curve(cook_lim, pn, h)
    bad <- is.na(c05) | is.nan(c05) | is.infinite(c05)
    for (i in which(bad)) {
      if (i > 1) {c05[i] <- c05[i - 1]}
    }
    
    fit <- loess(c05 ~ h, span = 0.2)
    c05 <- predict(fit)
    
    res_lim <- max(3.5,max(abs(df$stdt),na.rm = TRUE))
    g4<-plot_ly(data = df,  x = ~leverage,y = ~stdt,type = "scatter", color=~outlier,
            mode = "markers",colors= c( "GOOD" = "#A4D65E","outlier" = "#FF7F7F"),
            text = ~paste("idRow:", idRow,"<br>CookD:", round(cookD,4),"<br>Leverage:", round(leverage,4),"<br>Residuals:", round(stdt,4)),
            hoverinfo = "text",
            marker = list(size = pmax(6, sqrt(df$cookD) * 10)),
            showlegend=F
    ) %>%
      add_trace(x = h,y = c05,inherit=FALSE,type = "scatter",mode = "lines",
                name = paste0("CookD_thr = ",round(cook_lim,5)),line = list(color = "orange", dash = "dash",width = 2)
      ) %>%
      add_trace(x = h,y = -c05,inherit=FALSE,type = "scatter",mode = "lines",showlegend = FALSE,
                line = list(color = "orange",dash = "dash",width = 2)
      ) %>%
      add_trace(x = c(lev_lim, lev_lim),y = c(-3.5, 3.5),type = "scatter",mode = "lines",name = paste0("Leverage_thr= ",round(lev_lim,5)),
                line = list(color = "cyan",dash = "dash"),inherit = FALSE)%>%
      add_trace(x = c(0,lev_max),y = c(3.5, 3.5),type = "scatter",mode = "lines",name = "StdRes_thr = ±3.5",
                line = list(color = "blue",dash = "dot"),inherit = FALSE)%>%
      add_trace(x = c(0,lev_max),y = c(-3.5, -3.5),type = "scatter",mode = "lines",showlegend=F,
                line = list(color = "blue",dash = "dot"),inherit = FALSE)%>%
      layout(yaxis = list(range = c(-(res_lim+0.05), (res_lim+0.05))),
             xaxis = list(range = c(0, (lev_max+0.05))) )
    
    p<-subplot(g3, g4,g2, g1,nrows = 2,margin = 0.06) %>%  
      layout( height = 800,xaxis = list(title = list(text = "Residuals",font = list(size = 10)),tickfont = list(size = 8)),
        yaxis = list(title = list(text = "Frequency",font = list(size = 10)),tickfont = list(size = 8)),
        xaxis2 = list(title = list(text = "Leverage",font = list(size = 10)),tickfont = list(size = 8)),
        yaxis2 = list(title = list(text = "Residuals",font = list(size = 10)),tickfont = list(size = 8)),
        xaxis3 = list(title = list(text = "Theoretical Quantiles",font = list(size = 10)),tickfont = list(size = 8)),
        yaxis3 = list(title = list(text = "Sample Quantiles",font = list(size = 10)),tickfont = list(size = 8)),
        xaxis4 = list(title = list(text = "Fitted",font = list(size = 10)),tickfont = list(size = 8)),
        yaxis4 = list(title = list(text = "Raw Residuals",font = list(size = 10)),tickfont = list(size = 8)),
        annotations = list(
          list(text = "Histogram of StdResiduals",font = list(size = 13),x = 0.02,y = 1.025,xref = "paper",yref = "paper",showarrow = FALSE),
          list(text = "StdResiduals vs Leverage",font = list(size = 13),x = 0.71,y = 1.025,xref = "paper",yref = "paper",showarrow = FALSE),
          list(text = "QQ Plot",font = list(size = 13),x = 0.02,y = 0.45,xref = "paper",yref = "paper",showarrow = FALSE),
          list(text = "Raw Residuals vs Fitted",font = list(size = 13),x = 0.71,y = 0.45,xref = "paper",yref = "paper",showarrow = FALSE)
        )
      )
    
  return(p)
}
    