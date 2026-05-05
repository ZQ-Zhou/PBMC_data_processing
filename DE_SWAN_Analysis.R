age_cutoff_analysis <- function(data, min_age = 33, max_age = 70, max_samples = 10, step = 3) {
  results <- list()
  
  for (age_cutoff in seq(min_age, max_age, by = step)) {

    younger <- data %>% 
      filter(Age < age_cutoff) %>% 
      arrange(desc(Age))  
    
    older <- data %>% 
      filter(Age >= age_cutoff) %>% 
      arrange(Age)        
    
    n <- min(nrow(younger), nrow(older), max_samples)
    
    if (n < 2) {
      results[[as.character(age_cutoff)]] <- data.frame(
        Age_Cutoff = age_cutoff,
        Significant_Proteins = NA,
        N_Samples = n
      )
      next
    }

    younger <- younger %>% head(n)
    older <- older %>% head(n)

    p_values <- sapply(colnames(data)[-1], function(protein) {
      t.test(younger[[protein]], older[[protein]])$p.value
    })
    

    adj_p <- p_values
    

    sig_proteins <- sum(adj_p < 0.05, na.rm = TRUE)
    

    results[[as.character(age_cutoff)]] <- data.frame(
      Age_Cutoff = age_cutoff,
      Significant_Proteins = sig_proteins,
      N_Samples = n
    )
  }
  

  result_df <- bind_rows(results)
  
  return(result_df)
}
   
data = read.csv('D:/data/画图/PBMC-肺癌/data_merge_age_gender.csv',header = T)
data = data[,-2]
data_age = data %>% filter(label == 'H') %>% filter(Gender == '1') 
result2 <- age_cutoff_analysis(data_age[,-c(1,3)], min_age = 32, max_age = 70,step = 3)
data_age = data %>% filter(label == 'H') %>% filter(Gender == '0') 
result <- age_cutoff_analysis(data_age[,-c(1,3)], min_age = 32,max_age = 70,step = 3)

ggplot() +

     geom_line(
         data = result, 
         aes(x = Age_Cutoff, y = Significant_Proteins, color = "Female"), 
         linewidth = 2
       ) +

     geom_line(
         data = result2,
         aes(x = Age_Cutoff, y = Significant_Proteins, color = "Male"),
         linewidth = 2
       ) +

     geom_vline(
         xintercept = 50, 
         color = "#E15759", 
         linetype = "dashed", 
         linewidth = 1.5
       ) +

     scale_color_manual(
         name = NULL,
         values = c("Female" = "steelblue", "Male" = "darkorange"),
         guide = guide_legend(
             title.position = "top",
             title.hjust = 0.5,
             label.position = "right",
             keywidth = unit(1.5, "cm")
           )
       ) +

     annotate(
         "text", 
         x = 50, 
         y = max(c(result$Significant_Proteins, result2$Significant_Proteins), na.rm = TRUE) * 0.95,
         label = "Age 50", 
         color = "black", 
         size = 6, 
         vjust = 1.5,
         hjust = -0.1
       ) +

     labs(

           x = "Age",
         y = "Number of Significant Proteins (p < 0.05)",
         color = "Gender"
       ) +

     theme_bw() +
     theme(
         plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 24, face = "bold"),
         axis.text = element_text(size = 20),
         axis.text.y = element_text(size = 20, face = "bold"),
         axis.text.x = element_text(size = 20, face = "bold"),
         legend.text = element_text(size = 20),
         legend.position = c(0.35, 0.75),  
         legend.justification = c(1, 1),  
         legend.background = element_blank(),  
         legend.key = element_blank(),         
         legend.title = element_blank(), 
         panel.border = element_rect(fill = NA, linewidth = 1.5)
       ) +

     scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) 



