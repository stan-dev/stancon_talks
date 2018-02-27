library(Ryacas)
library(tidyverse)
library(stringr)

#functions for creating M-Splines as Ryacas expr. These are symbolic
#expressions that can be integrated symbolically by Ryacas, converted
#to strings, or converted to R functions.

M1 <- function(i,K,ts) {
  if(i < 4 | i > K + 4) {
    ret <- 0
  }
  else {
    #if its the same node avoid dividing by zero. just set to 0
    if(ts[[i+1]]-ts[[i]] == 0) ret <- 0
    else ret <- 1/(ts[[i+1]]-ts[[i]])
  }
  
  return(ret)
}

#recursion for creating higher order splines (greater than M=1)
M_ik <- function(i,k,ts,M1,M2,y) {
  if(ts[[i+k]]-ts[[i]] == 0) ret <- expression(0)
  else ret <- (k/(k-1))*((y-ts[[i]])*M1 + (ts[[i+k]]-y)*M2)/(ts[[i+k]]-ts[[i]])
  
  return(ret)
}

M2 <- function(i,ts,M1_expr,y) {
  
  A <- M_ik(i,2,ts,M1_expr[[i]],expression(0),y)
  B <- M_ik(i,2,ts,expression(0),M1_expr[[i+1]],y)
  
  return(list(A = A, B = B))
}

M3 <- function(i,ts,M2_expr,y) {
  A <- M_ik(i,3,ts,M2_expr[[i]]$A,expression(0),y)
  B <- M_ik(i,3,ts,M2_expr[[i]]$B,M2_expr[[i+1]]$A,y)
  C <- M_ik(i,3,ts,expression(0),M2_expr[[i+1]]$B,y)
  
  return(list(A=A, B=B,C=C))
}

M4 <- function(i,ts,M3_expr,y) {
  A <- M_ik(i,4,ts,M3_expr[[i]]$A,expression(0),y)
  B <- M_ik(i,4,ts,M3_expr[[i]]$B,M3_expr[[i+1]]$A,y)
  C <- M_ik(i,4,ts,M3_expr[[i]]$C,M3_expr[[i+1]]$B,y)
  D <- M_ik(i,4,ts,expression(0),M3_expr[[i+1]]$C,y)
  
  return(list(A=A, B=B,C=C,D=D))
}

#return Ryacas expr for all orders of M-Spline up to 4
#with K intermediary nodes. nodes param is a list in the
#with for example K=3 would be
#n <- list(t1=-20,t2=-20,t3=-20,t4=-20,t5=-10,t6=0,t7=10,t8=20,t9=20,t10=20,t11=20)
get_order_4_M_spline_expr <- function(K, nodes) {
  
  if(K != length(nodes) - 8) stop("nodes list does not have correct number of intermediary nodes, K")
  
  M <- 4
  y <- Sym("y")
  
  M1_expr <- map(1:(K+2*M-1), M1, K = K, ts = nodes)
  M2_expr <- map(1:(K+2*M-2),M2,ts=nodes,M1_expr = M1_expr, y = y)
  M3_expr <- map(1:(K+2*M-3),M3,ts=nodes,M2_expr = M2_expr, y = y)
  M4_expr <- map(1:(K+2*M-4),M4,ts=nodes,M3_expr = M3_expr, y = y)
  
  return(list(M1_expr = M1_expr, M2_expr = M2_expr, M3_expr = M3_expr, M4_expr = M4_expr))
}

#each function is represented as a list of expr because they are defined piecewise.
#first we must convert to an R function with a single argument which yacas_to_R_func does
#the create_MX_piecewise functions then take these list of R functions and puts them together
#to create a single R function with a single argument
yacas_to_R_func <- function(expr) {function(s) {expr %>% as.expression %>% eval(list(y=s))}}

#for the M1 splines we actually don't need to convert from yacas to R, because we don't save as expr
create_M1_piecewise <- function(f,i,e) {function(x) ifelse(between(x, e[[i]], e[[i+1]]), f, 0)}

create_M2_piecewise <- function(f_expr,i,e) {
  f <- map(f_expr, yacas_to_R_func)
  g <- function(x) { 
    if(x < e[[i]]) 0
    else if(between(x, e[[i]], e[[i+1]])) (f$A)(x)
    else if(between(x, e[[i+1]], e[[i+2]])) (f$B)(x)
    else 0
  }
  
  return(Vectorize(g))
}

create_M3_piecewise <- function(f_expr,i,e) {
  f <- map(f_expr, yacas_to_R_func)
  g <- function(x) { 
    if(x < e[[i]]) 0
    else if(between(x, e[[i]], e[[i+1]])) (f$A)(x)
    else if(between(x, e[[i+1]], e[[i+2]])) (f$B)(x)
    else if(between(x, e[[i+2]], e[[i+3]])) (f$C)(x)
    else 0
  }
  
  return(Vectorize(g))
}

create_M4_piecewise <- function(f_expr,i,e) {
  f <- map(f_expr, yacas_to_R_func)
  g <- function(x) { 
    if(x < e[[i]]) 0
    else if(between(x, e[[i]], e[[i+1]])) (f$A)(x)
    else if(between(x, e[[i+1]], e[[i+2]])) (f$B)(x)
    else if(between(x, e[[i+2]], e[[i+3]])) (f$C)(x)
    else if(between(x, e[[i+3]], e[[i+4]])) (f$D)(x)
    else 0
  }
  
  return(Vectorize(g))
}

#takes a list of functions, still in expr format (i.e a list of expr itself)
#and converts them to a list of usable R functions
list_of_expr_to_list_of_R_functions <- function(func_list_expr, piecewise_func_creator, nodes) {
  map2(func_list_expr, 1:length(func_list_expr), piecewise_func_creator, e = nodes)
}

#easy outward facing function to conver output of get_order_4_M_spline_expr
#to list of list of usable R functions
expr_to_R_functions <- function(expr_list, nodes) {
  map2(expr_list,
       list(create_M1_piecewise,create_M2_piecewise,create_M3_piecewise,create_M4_piecewise),
       list_of_expr_to_list_of_R_functions,
       nodes = nodes)
}

#function to take list of R functions and plot. useful for visualizing splines
plot_function_list <- function(fl, nodes) {
  node_df <- tibble(node = names(nodes), location = unlist(nodes))
  colors <- c("#F8766D", "#00BFC4", "#B79F00", "#619CFF", "#00BA38", "#F564E3")
  p <- ggplot(data.frame(x = c(min(unlist(nodes)), max(unlist(nodes)))), aes(x))
  for(i in 1:length(fl)) p <- p + stat_function(fun = fl[[i]], geom = "line", color = colors[[(i %% 6) + 1]])
  p + geom_label(aes(x = location, y = 0, label = node), data = node_df)
}

#convert M4_expr to I-Spline
symbolic_Mspline_to_R_Ispline <- function(i, M4_expr, nodes) {
  e <- nodes
  y <- Sym("y")
  expr_list <- map(M4_expr[[i]], function(expr) as.expression(Integrate(yacas(expr),y)))
  
  f <- function(s) {
    ret <- 0
    
    if(s < e[[i]]) return(0)
    
    if(between(s,e[[i]],e[[i+1]])) return(ret + eval(expr_list[[1]], list(y=s)) - eval(expr_list[[1]], list(y=e[[i]])))
    else ret <- ret + eval(expr_list[[1]], list(y=e[[i+1]])) - eval(expr_list[[1]], list(y=e[[i]]))
    
    if(between(s,e[[i+1]],e[[i+2]])) return(ret + eval(expr_list[[2]], list(y=s)) - eval(expr_list[[2]], list(y=e[[i+1]])))
    else ret <- ret + eval(expr_list[[2]], list(y=e[[i+2]])) - eval(expr_list[[2]], list(y=e[[i+1]]))
    
    if(between(s,e[[i+2]],e[[i+3]])) return(ret + eval(expr_list[[3]], list(y=s)) - eval(expr_list[[3]], list(y=e[[i+2]])))
    else ret <- ret + eval(expr_list[[3]], list(y=e[[i+3]])) - eval(expr_list[[3]], list(y=e[[i+2]]))
    
    if(between(s,e[[i+3]],e[[i+4]])) return(ret + eval(expr_list[[4]], list(y=s)) - eval(expr_list[[4]], list(y=e[[i+3]])))
    else ret <- ret + eval(expr_list[[4]], list(y=e[[i+4]])) - eval(expr_list[[4]], list(y=e[[i+3]]))
    
    if(s > e[[i+4]]) ret <- 1
  
    return(ret)
  }
  
  return(Vectorize(f))
}

#takes list of M4 functions as expr and convert to list of R functions
M4_expr_to_R_I_spline <- function(M4_expr, nodes) {
  map(1:length(M4_expr), symbolic_Mspline_to_R_Ispline, M4_expr = M4_expr, nodes = nodes)
}

#create a function from a linear combination of the I-Splines
convex_combination <- function(simplex, I_splines) {
  #functions <- map(1:length(simplex), function(n) symbolic_Mspline_to_R_Ispline(n,M4_expr))
  f <- function(x) reduce(map2(I_splines, simplex, function(f,s) s*f(x)),sum)
  
  return(Vectorize(f))
}

#take expression and translate to Stan code
"%+%" = function(x,y) {
  if(is.character(x) || is.character(y)) {
    return(paste(x , y, sep=""))
  } else {
    .Primitive("+")(x,y)
  }
}

#take the M4 expressions, integrates and returns Stan code of the I-Splines
expr_to_stan_str <- function(i, M4_expr, nodes) {
  e <- nodes
  y <- Sym("y")
  expr_list <- map(M4_expr[[i]], function(expr) as.expression(Integrate(yacas(expr),y)))
  
  expr_functions <- map(expr_list, as.expression) %>% map(function(expr) {function(x) eval(expr, list(y = x))})
  
  expr_strings <- expr_list %>%
    map(as.character) %>%
    str_replace("Integrate\\( y \\) ", "") %>%
    str_replace_all("y", "x")

  str <- "real I" %+% as.character(i) %+% "(real x) { \n"
  str <- str %+% "\t real ret = 0;\n"
  str <- str %+% "\t if(x < " %+% as.character(e[[i]]) %+% ") return(0);\n\n"
  
  str <- str %+% "\t if(x >= " %+% as.character(e[[i]]) %+% " && x < " %+% as.character(e[[i+1]]) %+% ") return(ret + " %+% expr_strings[[1]] %+% " - " %+% as.character(expr_functions[[1]](e[[i]])) %+% ");\n"
  str <- str %+% "\t else ret = ret + " %+% as.character(expr_functions[[1]](e[[i+1]]) - expr_functions[[1]](e[[i]])) %+% ";\n\n"
  
  str <- str %+% "\t if(x >= " %+% as.character(e[[i+1]]) %+% " && x < " %+% as.character(e[[i+2]]) %+% ") return(ret + " %+% expr_strings[[2]] %+% " - " %+% as.character(expr_functions[[2]](e[[i+1]])) %+% ");\n"
  str <- str %+% "\t else ret = ret + " %+% as.character(expr_functions[[2]](e[[i+2]]) - expr_functions[[2]](e[[i+1]])) %+% ";\n\n"
  
  str <- str %+% "\t if(x >= " %+% as.character(e[[i+2]]) %+% " && x < " %+% as.character(e[[i+3]]) %+% ") return(ret + " %+% expr_strings[[3]] %+% " - " %+% as.character(expr_functions[[3]](e[[i+2]])) %+% ");\n"
  str <- str %+% "\t else ret = ret + " %+% as.character(expr_functions[[3]](e[[i+3]]) - expr_functions[[3]](e[[i+2]])) %+% ";\n\n"
  
  str <- str %+% "\t if(x >= " %+% as.character(e[[i+3]]) %+% " && x < " %+% as.character(e[[i+4]]) %+% ") return(ret + " %+% expr_strings[[4]] %+% " - " %+% as.character(expr_functions[[4]](e[[i+3]])) %+% ");\n"
  str <- str %+% "\t else ret = ret + " %+% as.character(expr_functions[[4]](e[[i+4]]) - expr_functions[[4]](e[[i+3]])) %+% ";\n\n"
  
  str <- str %+% "\t if(x >= " %+% as.character(e[[i+4]]) %+% ") return(1);\n\n"
  
  str <- str %+% "\t return(ret);\n"
  str <- str %+% "}"
  
  return(str)
}

#map(1:length(M4_expr), expr_to_stan_str, M4_expr = M4_expr) %>% paste(collapse = "\n\n") %>% cat

#create convex combination function in Stan
linear_combination_stan <- function() {
  str <- "real[] ispline(real[] x, matrix c) {\n\n"
  str <- str %+% "\t int N;\n"
  str <- str %+% "\t real ret[dims(x)[1]];\n"
  str <- str %+% "\t N = dims(x)[1];\n\n"
  str <- str %+% "\t for(n in 1:N) ret[n] = "
  for(i in 1:(length(M4_expr)-1)) str <- str %+% "c[n][" %+% as.character(i) %+% "]*I" %+% as.character(i) %+% "(x[n]) + "
  str <- str %+% "c[n][" %+% as.character(length(M4_expr)) %+% "]*I" %+% as.character(length(M4_expr)) %+% "(x[n]);\n\n"
  str <- str %+% "\t return(ret);\n}"

  return(str)
}
