library(DT)
tab1 <- read.csv('reports/genomestats.csv')
tab2 <- DT::datatable(tab1, class = 'order-column stripe', options = list(pageLength = 50), rownames = F)
saveWidget(tab2, 'reports/genomestats.html', selfcontained = TRUE)

tab3 <- read.table('reports/checkm.summary', header = TRUE, sep = "\t", comment.char = "")
tab4 <- DT::datatable(tab3, class = 'order-column stripe', options = list(pageLength = 50), rownames = F)
saveWidget(tab4, 'reports/checkm.summary.html', selfcontained = TRUE)
