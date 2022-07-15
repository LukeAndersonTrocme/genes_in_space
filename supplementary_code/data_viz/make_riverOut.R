make_riverOut <- function(df) {
  # filial links
  nodes.x <-
    df %>%
    distinct(gen, yax, rgb, labs) %>% 
    dplyr::rename(yax = yax, gen = gen, rgb = rgb, labels = labs) 
  # parental links
  nodes.y <-
    df %>% 
    distinct(gen_p, yax_mp, rgb_mp) %>% 
    dplyr::rename(yax = yax_mp, gen = gen_p,  rgb = rgb_mp) %>%
    mutate(labels = "")
  
  nodes <- 
    rbind(nodes.x, nodes.y) %>% unique() %>% 
    dplyr::mutate(ID = paste(gen,yax, sep = "_"),
                  x = gen#,
                  #y = dense_rank(yax)
    )
  
  unique_y <- 
    nodes %>% distinct(yax) %>% arrange(yax) %>% 
    dplyr::mutate(y = row_number())
  
  nodes <- 
    inner_join(nodes, unique_y) %>%
    dplyr::select(x, y, ID, labels, rgb)%>% 
    as.data.frame(stringsAsFacors=FALSE)
  
  node_styles <- 
    lapply(nodes$ID, function(n) 
    {list(col = nodes[which(nodes$ID == n),]$rgb, 
          nodestyle = "point", 
          #textcol="black",
          #cex = 5,
          #pos="left",
          srt = 90 #70
    )})
  names(node_styles) = nodes$ID
  #re-scale line width
  scale_range <- function(x){(x-min(x))/(max(x)-min(x))/5}
  #links are per migration event
  edges <-
    df %>%
    mutate(N1 = paste( gen_p, yax_mp, sep = "_"),
           N2 = paste(gen, yax, sep = "_"),
           Value = scale_range(total_contribution)) %>% ######### HERE #########
  dplyr::select(N1, N2, Value) %>%
    as.data.frame() %>% filter(N2 %in% nodes$ID)
  
  riverOut <-
    makeRiver(nodes = nodes,
              edges = edges,
              node_styles = node_styles)
}