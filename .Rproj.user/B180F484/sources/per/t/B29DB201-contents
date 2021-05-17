#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Created by Joonho Gong
# Email : jgong9@ncsu.edu


## For Shiny and plots
library(shiny)
library(shinythemes)
library(DT)
library(latex2exp)
library(shinydashboard)
library(markdown)


## For model fitting
library(readxl)
library(minpack.lm)
library(nlstools)


##### Read in the original data file
isotherm.meta<-read_excel("data/Glyphosate_Paraquat_master.xlsx",sheet="isotherm.meta")
isotherm.data <- read_excel("data/Glyphosate_Paraquat_master.xlsx",sheet="isotherm.data",skip=10, .name_repair = "minimal")
data_set <- isotherm.data
data_meta <- data.frame(isotherm.meta)

n_datasets <- ncol(isotherm.data) / 2

x_matrix <- matrix(0, n_datasets, 101)
n_vec <- vector("numeric", length=n_datasets)


##### Functions required
## 1. Langmuir
negLogLikelihood <- function(par, data_x, data_y) {

  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(data_y, mean = (par[1] * data_x * exp(par[2])) / (1 +  exp(par[2]) *data_x) , sd = par[3])
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle <- function(Q, K_new, C, sigma, cov_mat){
  der_vec <- rep(NA, 3)
  der_vec[1] <- (exp(K_new) * C) / (1 + exp(K_new) * C)
  der_vec[2] <- (Q  * exp(K_new) * C) / (1 + exp(K_new) * C)^(2) 
  der_vec[3] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec  ) 
  
}


## 2. Freundlich
negLogLikelihood_fre <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/par[2]) ), sd = par[3])
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

negLogLikelihood_fre_trans <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(data_y, mean = (exp(par[1]) * data_x^(1/par[2]) ), sd = par[3])
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle_fre <- function(K_F, n, C, sigma, cov_mat){
  der_vec <- rep(NA, 3)
  der_vec[1] <- C^(1/n)
  der_vec[2] <- K_F * log(C) * C^(1/n) * (-1/n^(2))
  der_vec[3] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec ) 
  
}

se_for_CI_mle_fre_trans <- function(K_F_prime, n, C, sigma, cov_mat){
  der_vec <- rep(NA, 3)
  der_vec[1] <- exp(K_F_prime) * C^(1/n)
  der_vec[2] <- exp(K_F_prime) * log(C) * C^(1/n) * (-1/n^(2))
  der_vec[3] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec ) 
  
}
## Hessian function
# hes_fnt <- function(par, data_x, data_y){
#     result_mat <- matrix(0, length(par), length(par))
#     desing_mat <- cbind( rep(1, length(data_x)), data_x)
#     result_mat[1:(length(par)-1), 1:(length(par)-1)] <- - 1/par[3]^(2) *  t(desing_mat) %*% desing_mat
#     -1/ (par[3])^(2) * data_x
#     
#     side <- -1/ (par[3])^(4) * t(desing_mat) %*% (data_y - desing_mat %*% c(par[1], par[2]))
#     result_mat[1:(length(par)-1), length(par)] <- side
#     result_mat[length(par), 1:(length(par)-1)] <- t(side)
#     result_mat[length(par),length(par)] <- length(data_x)/ (2 * par[3]^(4)) -1/par[3]^(6) * t((data_y - desing_mat %*% c(par[1], par[2]))) %*% (data_y - desing_mat %*% c(par[1], par[2]))
#             ## the function to minimize is NEGATIVE log likelihood, so need -1
#     return( -result_mat)
# }

##############
## 3. Sips ###
############## 
## Original
# negLogLikelihood_sips <- function(par, data_x, data_y) {
#   
#   
#   if(par[3] < 0){
#     output <- 10000000
#   } else {
#     # likelihoods <- dnorm(data_y, mean = (par[1]  ) / (1 +  ( exp(par[2]) * data_x^(1/ par[3])^(-1)  ) ) , sd = par[4])
#     likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/par[3] ) * exp(par[2])) / (1 +  exp(par[2]) * data_x^(1/ par[3]  ) ) , sd = par[4])
#     
#     log.likelihoods <- log(likelihoods)
#     output <- -sum(log.likelihoods)
#     return(output)
#   }
# } # end of negLogLikelihood
# 
# se_for_CI_mle_sips <- function(Q, K_new, n, C, sigma, cov_mat){
#   der_vec <- rep(NA, 4)
#   der_vec[1] <- (exp(K_new) * C^(1/ n   ) ) / (1 + exp(K_new) * C^(1/ n  ) )
#   der_vec[2] <- (Q  * exp(K_new) * C^(1/ n  ) ) / (1 + exp(K_new) * C^(1/ n  ) )^(2) 
#   der_vec[3] <- Q  * exp(K_new) * -(C^(-1/ n ) + exp(K_new))^(-2) * C^(-1/n  ) * log(C) * (n )^(-2)
#   der_vec[4] <- 0
#   
#   se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec  ) 
#   
# }

negLogLikelihood_sips <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
    # likelihoods <- dnorm(data_y, mean = (par[1]  ) / (1 +  ( exp(par[2]) * data_x^(1/ par[3])^(-1)  ) ) , sd = par[4])
    ## short form
    # likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(par[3] ) * exp(par[2])) / (1 +  exp(par[2]) * data_x^(par[3]  ) ) , sd = par[4])
    ## revse 1/n
    # likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/par[3] ) * par[2]) / (1 +  par[2] * data_x^( 1/par[3]  ) ) , sd = par[4])
    ## no trans at all
    likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/par[3] ) * exp(par[2]) ) / (1 +  exp(par[2]) * data_x^( 1/par[3]  ) ) , sd = par[4])
    
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle_sips <- function(Q, K_new, n, C, sigma, cov_mat){
  der_vec <- rep(NA, 4)
  der_vec[1] <- (exp(K_new) * C^(1/n   ) ) / (1 + exp(K_new) * C^(1/n  ) )
  der_vec[2] <- (Q  * exp(K_new) * C^( 1/n  ) ) / (1 + exp(K_new) * C^(1/n  ) )^(2)
  der_vec[3] <- Q  * exp(K_new) * -(C^(-1/n ) + exp(K_new))^(-2) * C^(-1/n  ) * log(C) * (n)^(-2)
  der_vec[4] <- 0
  
  # der_vec <- rep(NA, 4)
  # der_vec[1] <- (exp(K_new) * C^(n   ) ) / (1 + exp(K_new) * C^(n  ) )
  # der_vec[2] <- (Q  * exp(K_new) * C^( n  ) ) / (1 + exp(K_new) * C^(n  ) )^(2) 
  # der_vec[3] <- Q  * exp(K_new) * -(C^(-n ) + exp(K_new))^(-2) * C^(-n  ) * log(C) * -1
  # der_vec[4] <- 0
  ## reverse 1/n
  
  # der_vec <- rep(NA, 4)
  # der_vec[1] <- (K_new * C^(1/n   ) ) / (1 + K_new * C^(1/n  ) )
  # der_vec[2] <- (Q  * K_new * C^( 1/n  ) ) / (1 + K_new * C^(1/n  ) )^(2) 
  # der_vec[3] <- Q  * K_new * -(C^(-1/n ) + K_new)^(-2) * C^(-1/n  ) * log(C) * (n)^(-2)
  # der_vec[4] <- 0
  # ## no trans at all
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec  ) 
  
}



negLogLikelihood_sips_trans <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
     likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/ (1+exp(par[3])) ) * exp(par[2])) / (1 +  exp(par[2]) * data_x^(1/ (1+exp(par[3]) )  ) ) , sd = par[4])
    ## original
    # likelihoods <- dnorm(data_y, mean = (par[1]  ) / (1 + ( exp(par[2]) * data_x^(1/ (1+exp(par[3]) )  ) )^(-1) ) , sd = par[4])
    ## short with -1
    # likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/ (exp(par[3])) ) * exp(par[2])) / (1 +  exp(par[2]) * data_x^(1/ (exp(par[3]) )  ) ) , sd = par[4])
    ## without 1+
    
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle_sips_trans <- function(Q, K_new, n_new, C, sigma, cov_mat){
  der_vec <- rep(NA, 4)
  der_vec[1] <- (exp(K_new) * C^(1/ (1+exp(n_new) )   ) ) / (1 + exp(K_new) * C^(1/ (1+exp(n_new) )  ) )
  der_vec[2] <- (Q  * exp(K_new) * C^(1/ (1+exp(n_new) )  ) ) / (1 + exp(K_new) * C^(1/ (1+exp(n_new) )  ) )^(2)
  der_vec[3] <- Q  * exp(K_new) * -(C^(-1/ (1+exp(n_new) ) ) + exp(K_new))^(-2) * C^(-1/ (1+exp(n_new) )  ) * log(C) * (1+exp(n_new) )^(-2) * exp(n_new)
  # trans original
  
  # der_vec[1] <- (exp(K_new) * C^(1/ (exp(n_new) )   ) ) / (1 + exp(K_new) * C^(1/ (exp(n_new) )  ) )
  # der_vec[2] <- (Q  * exp(K_new) * C^(1/ (exp(n_new) )  ) ) / (1 + exp(K_new) * C^(1/ (exp(n_new) )  ) )^(2) 
  # der_vec[3] <- Q  * exp(K_new) * -(C^(-1/ (exp(n_new) ) ) + exp(K_new))^(-2) * C^(-1/ (exp(n_new) )  ) * log(C) * C^(-n_new)
  # ## trans without 1+
  
  der_vec[4] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec  ) 
}

negLogLikelihood_sips_linear <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(1/data_y, mean = par[1] + par[1] * par[2] * data_x^( 1/ par[3] ) , sd = par[4])
    # likelihoods <- dnorm(1/data_y, mean = par[1] + par[2] * data_x^( 1/ (1+exp(par[3]) ) ) , sd = par[4])
    
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle_sips_linear <- function(Q, K, n, C, sigma, cov_mat){
  der_vec <- rep(NA, 4)
  der_vec[1] <- 1 + K * C^(1/n)
  der_vec[2] <- Q * C^(1/n)
  der_vec[3] <- Q * K * (-n) * C^(-n-1)
  der_vec[4] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec  ) 
  
}


sidebar <- dashboardSidebar(
  sidebarMenu(id="tabs",
              menuItem("Initial values", tabName="initial", icon=icon("file-text-o"), selected=T),
              menuItem("Confidence envelopes", tabName = "conf", icon=icon("line-chart"), selected=F),
              menuItem("AIC", tabName = "aic", icon=icon("pencil"), selected=F)
            #  menuItem("Question", tabName = "ques", icon = icon("question"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "initial",
              br(),

            
              fluidRow(
                column(
                  p(
                    strong("1. Langmuir model"),
                    br(),
                    br(),
                       "In order to get the initial values for \\(Q_{max}\\), \\(K_{d}\\), and \\(\\sigma^{2}\\), we transform them as 
                       \\(Q^{\\star}_{max} = Q_{max}^{-1}\\), \\(K^{\\star}_{d} = K_{d}^{-1}\\). Let \\(Y^{\\star}=q^{-1}\\) and \\(C^{\\star}_{W}= C^{-1}_{W}\\). Then, ",
                       br(),
                       br(),
                       " \\( q = \\frac{Q_{max} K_{d} C_{W}}{1 + K_{d} C_{W}} \\Longleftrightarrow q = \\frac{Q_{max} C_{W}}{K_{d}^{\\star} +  C_{W}}\\)",
                       br(),
                      "\\(\\Longrightarrow \\frac{1-q}{q} = \\frac{K_{d}^{\\star} + C_{W} - Q_{max}C_{W}}{Q_{max}C_{W}}\\)",
                      br(),
                      "\\(\\Longleftrightarrow \\frac{1-q}{q} = Q_{max}^{\\star} + K_{d}^{\\star} Q_{max}^{\\star} C_{W}^{\\star} - 1\\)",
                      br(),
                      "\\(\\Longrightarrow Y^{\\star} =  Q_{max}^{\\star} + K_{d}^{\\star} Q_{max}^{\\star} C_{W}^{\\star}\\).",
                      br(),
                    br(),
                      "Since we now have a linear model with an intercept \\(Q_{max}^{\\star}\\) and a slope \\(K_{d}^{\\star}Q_{max}^{\\star}\\), we can apply the ordinary least square (OLS) method to estimate both of them and get the initial values. 
                      Specifically, we set \\(\\hat{Q}_{max} = \\frac{1}{\\hat{Q}^{\\star}_{max}}\\) and \\(\\hat{K}_{d}=\\frac{\\hat{Q}^{\\star}_{max}}{\\widehat{K_{d}^{\\star} Q_{max}^{\\star}}}\\). For \\(\\sigma^{2}\\), we adopt the sum of squared errors \\(\\hat{\\sigma}^{2} = \\sum^{n}_{i=1} ( y^{\\star}_{i} - \\hat{y}^{\\star}_{i}  ). \\)", 
                      style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"
                    ),
                  br(),
                  p(
                    strong("2. Freundlich model"),
                    br(),
                    br(),
                    "Similarly, we reparameterize the Freundlich model with \\(q^{\\star}=\\log{q}\\), \\(C^{\\star}_{W}=\\log{C_{W}}\\), \\(K^{\\star}_{F} = \\log{K_{F}}\\), and \\(n^{\\star} = \\frac{1}{n}\\) . Then, ",
                    br(),
                    br(),
                    "\\( q^{\\star} = K^{\\star}_{F} + n^{\\star} C^{\\star}_{W}\\).",
                    br(),
                    br(),
                    "Since we now have a linear model with an intercept \\( K^{\\star}_{F}\\) and a slope \\(n^{\\star}\\), we can apply OLS again to estimate both of them and get the initial values. 
                      In particular, we set \\(\\hat{K}_{F} = exp(\\hat{K}^{\\star}_{F}) \\) and \\(\\hat{n}= \\frac{1}{\\hat{n}^{\\star} }\\). For \\(\\sigma^{2}\\), we adopt the sum of squared errors \\(\\hat{\\sigma}^{2} = \\sum^{n}_{i=1} ( q^{\\star}_{i} - \\hat{q}^{\\star}_{i}  ). \\)", 
                    style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"
                ),
                  width=11),
            ),
    ), # end of tabItme 1
    tabItem(tabName = "conf",
            br(),
            
            fluidRow(
              column(
                p(
                  strong("Delta method"),
                  br(),
                  br(),
                  "Let \\( \\mathbf{\\theta}_{0} = [Q_{max}, \\: K^{\\prime}_{d}, \\: \\sigma^{2} ]\\) and denote its maximum likelihood estimator (MLE) by \\( \\hat{\\mathbf{\\theta}}_{MLE}\\). Since the aymptotic distribution of MLE is multivariate normal, we know ",
                  "$$ \\sqrt{n}( \\hat{\\mathbf{\\theta}}_{MLE} - \\mathbf{\\theta}_{0} ) \\; \\overset{\\cdot}{\\sim}\\; N( \\mathbf{0}, \\mathbf{I}^{-1}(\\mathbf{\\theta}_{0}) ))$$",
                  "where \\(\\mathbf{I}\\) is the Fisher information matrix. Next, define a continuous function \\( g(\\mathbf{\\theta}_{0})=\\frac{Q_{max}exp(K^{\\prime}_{d}) C_{W} }{1+exp(K^{\\prime}_{d})C_{W} }\\). Since its partial derivaties exist, the delta method is applicable and 
                  we can get the \\((1-\\alpha)\\)% confidence envelops as below.",
                  "$$ g(\\hat{\\mathbf{\\theta}}_{MLE}) \\pm z_{\\frac{\\alpha}{2}} \\times \\sqrt{ \\mathbf{g}^{\\prime T}(\\hat{\\mathbf{\\theta}}_{MLE}) \\mathbf{H}^{-1}(\\hat{\\mathbf{\\theta}}_{MLE}) \\mathbf{g}^{\\prime}(\\hat{\\mathbf{\\theta}}_{MLE})  }, $$",
                  "where \\(\\mathbf{g}^{\\prime}\\) is derivative of the function \\(g\\), \\(\\mathbf{H}\\) is numerical hessian matrix, and \\(\\alpha\\) is a given significance level with corresponding critical value \\(z\\). ",
                  br(),
                  br(),
                  "Similarly, the confidence envelops for the Freundlich model can be obtained through the delta method with a continuous function \\(h(\\mathbf{\\theta})=exp(K_{F}) C^{\\frac{1}{n}}_{W}\\).",
                  style="text-align:justify;color:black;background-color:#D5F5E3;padding:15px;border-radius:10px"
                ),
                width=11),
            )

    ), # end of tabItem 2
    tabItem(tabName = "aic",
            fluidRow(
              column(
                br(),

                
                p(
                  strong("Akaike information criterion (AIC)"),
                  br(),
                  br(),
                  "Recall the models with Gaussian random error. Let k be the number of parameters in each model and \\(L\\)\ the likelihood function for the model. Then the AIC value of the model is defined as ",
                  "$$ \\text{AIC} = 2k - \\ln{L(\\hat{\\mathbf{\\theta}}_{MLE} )}.$$",
                  style="text-align:justify;color:black;background-color:#D8E8F1;padding:15px;border-radius:10px"
                ),
                width=11),
            )
    )
    
  ), # end of tabItems
) # end of dashboardBody




ui <- tagList(
  tags$style("@import url(https://use.fontawesome.com/releases/v5.6.0/css/all.css);"),
  
 # tags$head(HTML("<title>Toxicant Analysis</title>")),
  tags$head(HTML("<title>Toxicology</title> <link rel='icon' type=image/gif/png href='warning.png'>")),
  
  navbarPage(
        title = span("Toxin Analysis", style="font-weight: bold; font-weight: 900;"), 
        theme = shinytheme("flatly"),
       # shinythemes::themeSelector(),
       tabPanel("Model",
                style='padding-left:10px; padding-right:25px; padding-top:10px; padding-bottom:5px',
                
                fluidRow(
                  column(
                    br(),
                    div( img(src="acs-logo.svg",width="169px",height="54px"),  align= "center"),
                    p("For more information, please see ",em("Wang et al. (2019)"),
                      br(),
                      a(href="https://pubs.acs.org/doi/10.1021/acsomega.9b02051", "Here",target="_blank",style = "color: blue;"),style="text-align:center;color:black"),
                    br(),
                    br(),
                    br(),
                    div(img(height = 50, width = 140, src = "tam_logo2.png"),
                        strong("X"),
                        img(height = 50, width = 140, src = "ncstate_logo.png"),
                        align= "center"
                    ),
                    br(),
                    div("This is a product of   ", 
                        
                        a(actionButton(inputId = "email1", label = "Joonho Gong", 
                                       icon = icon("envelope", lib = "font-awesome"),
                                       style="color: #0303FC; background-color: #EBEDEF; border-color: #DCDCDC"
                        ),
                        href="mailto:jgong9@ncsu.edu"),
                        align= "center"
                    ),
                    width=3),
                  
                  
                  column(
                    h4(strong("Isotherm Models")),
                    p( strong("1. Langmuir model"),
                       withMathJax("$$ q = \\frac{Q_{max} K_{d} C_{W}}{1+K_{d} C_{W} },$$"),
                       "where \\(q\\) is toxin absorbed (mol/kg), \\(Q_{max}\\) is the maximum capacity (mol/kg), \\(K_{d}\\) is a distribution constant, and \\(C_{W}\\) is the toxin equilibirum conncentration.",style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                    br(),
                    
                    p( strong("2. Freundlich model"),
                       withMathJax("$$ q = K_{F} \\cdot C_{W}^{\\frac{1}{n}} $$"),
                       #   withMathJax("$$ q = \\frac{Q_{max} K_{d} C}{1+K_{d} C },$$"),
                       
                       "where \\(q\\) is toxin absorbed (mol/kg), \\(K_{F}\\) is the Freundlich constant, \\(\\frac{1}{n}\\) is the measure of intensity, and \\(C_{W}\\) is the toxin equilibirum conncentration.",
                       style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                    br(),
                    # p(strong("3. Sips model"),
                    #   withMathJax("$$ q = \\frac{Q_{max} K_{d} C^{\\frac{1}{n}} }{1+K_{d} C^{\\frac{1}{n}} },$$
                    #                \n where \\( 0 < \\frac{1}{n}\\le 1 \\)."),
                    #   '\n',
                    #   "Details about the model including what the parameters indicate.",
                    #   style="text-align:justify;color:black;background-color:#DFFEF7;padding:15px;border-radius:10px"),
                    width=8)
                  
                )
       ), ## end of tab 1
        tabPanel("Plot",
                 sidebarLayout(
                     sidebarPanel(
                         selectizeInput('main', 'Toxicant', c("Glyphosate", "Paraquat")
                                        ),
                      
                         conditionalPanel(
                             condition = "input.main=='Glyphosate'",
                             sliderInput("k",
                                         "Data ID of Glyphosate:",
                                         min = 1,
                                         max = 8,
                                         value = 1,
                                         step=1,
                                         round =T),
                         ),
                         conditionalPanel(
                             condition = "input.main=='Paraquat'",
                             sliderInput("k2",
                                         "Data ID of Paraquat:",
                                         min = 9,
                                         max = 12,
                                         value = 9,
                                         step=1,
                                         round=T),
                         ),
                         
                         selectizeInput('Rfunction', 'Choose your R optimization function', c("nlm", "optim")
                                        ),
                         
                         conditionalPanel(
                           condition = "input.Rfunction=='optim'",
                           selectizeInput('method', 'Method of minimization', c( "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"  ))
                         ),

                        selectizeInput("conf_level", "Confidence level for MLE band", c("99%", "95%", "90%", "80%")),
                        
                         br(),
                         br(),
                         br(),
                         br(),
                         br()
                         
                        #  div(img(height = 30, width = 100, src = "tam_logo2.png"),
                        #      "X",
                        #      img(height = 30, width = 100, src = "ncstate_logo.png"),
                        #      align= "center"
                        #  ),
                        # br(),
                        #  div("This is a product of ", 
                        # 
                        #      a(actionButton(inputId = "email1", label = "Joonho Gong", 
                        #                     icon = icon("envelope", lib = "font-awesome"),
                        #                     style="color: #0303FC; background-color: #EBEDEF; border-color: #DCDCDC"
                        #                     ),
                        #        href="mailto:jgong9@ncsu.edu"),
                        #      align= "center"
                        #  )
                     ),
                     mainPanel(
                                h3(strong("1. Langmuir Model")),
                                h4("1-1. Non-transformed OLS vs Transformed MLE"),
                                plotOutput("trajectoryPlot"),
                                withMathJax("For MLE, we assume \n
                                     $$ q = \\frac{Q_{max} \\cdot \\exp{K_{d}^{\\prime}} \\cdot C_{W}}{(1+ \\exp{K_{d}^{\\prime}} \\cdot C_{W})} + \\varepsilon, $$
                                     \n where \\(K_{d} = \\exp(K_{d}^{\\prime})\\) and \\( \\varepsilon \\sim N(0, \\sigma^{2}) \\).
                                     \n"),
                                br(),
                                br(),
                                h4("1-2. MLE summary of Langmuir"),
                                tableOutput("mle_table"),
                        
                                h4("1-3. AIC of Langmuir"),
                                verbatimTextOutput("aic_lan"),
                                # h4("1-3. Numerical Hessian of MLE and its invserse of Langmuir"),
                                # fluidRow(
                                #   column(6,
                                #          "a. Numerical Hessian of MLE",
                                #          tableOutput("hessian_num_original")
                                #         ),
                                #   column(6,
                                #          "b. Inverse of Hessian of MLE",
                                #          tableOutput("hessian_num_original_inverse")
                                #          )
                                #  ), # end of fluidRow
                                br(),
                                br(),
                                h3(strong("2. Freundlich vs Langmuir")),
                                h4("2-1. Trajectory plot"),
                                plotOutput("plot_fre"),
                                withMathJax("For MLE, we assume \n
                                     $$ q = \\exp{K^{\\prime}_{F}} \\cdot C_{W}^{\\frac{1}{n}}+ \\varepsilon, $$
                                     \n where \\(K_{F} = \\exp{K^{\\prime}_{F}}\\) and \\( \\varepsilon \\sim N(0, \\sigma^{2}) \\).
                                     \n"),
                                br(),
                                br(),
                                h4("2-2. MLE summary of Freundlich"),
                                tableOutput("mle_table_fre"),
                                h4("2-3. AIC of Freundlich"),
                                verbatimTextOutput("aic_fre"),
                                # h4("2-3. Numerical Hessian of MLE and its invserse of Freundlich"),
                                # fluidRow(
                                #   column(6,
                                #          "a. Numerical Hessian of MLE",
                                #          tableOutput("hessian_num_original_fre")
                                #   ),
                                #   column(6,
                                #          "b. Inverse of Hessian of MLE",
                                #          tableOutput("hessian_num_original_inverse_fre")
                                #   )
                                # ), # end of fluidRow
                                br(),
                                br(),
                                # h4("3. Sips vs Langmuir"),
                                # h4("3-1. Non-transformed OLS vs Transformed MLE"),
                                # plotOutput("plot_sips"),
                                # withMathJax("For MLE, we assume \n
                                #      $$ q = \\frac{Q_{max} \\cdot \\exp{K^{\\prime}} \\cdot C^{\\frac{1}{n}} }{(1+ \\exp{K^{\\prime}} \\cdot C^{\\frac{1}{n}} )} + \\varepsilon, $$
                                #      \n where \\(K = \\exp(K^{\\prime})\\) and \\( \\varepsilon \\sim N(0, \\sigma^{2}) \\).
                                #      \n"),
                                # br(),
                                # br(),
                                # h4("3-2. MLE summary of Sips"),
                                # tableOutput("mle_table_sips"),
                                # br(),
                                # h4("3-3. Inverse of numerical Hessian of Sips"),
                                # tableOutput("hessian_num_original_inverse_sips"),
                                # br(),
                                # h4("3-4. Transformed Sips"),
                                # plotOutput("plot_sips_trans"),
                                # withMathJax("For MLE, we assume \n
                                #      $$ q = \\frac{Q_{max} \\cdot \\exp{K^{\\prime}} \\cdot C^{\\frac{1}{1+\\exp{n^{\\prime} } }} }{(1+ \\exp{K^{\\prime}} \\cdot C^{\\frac{1}{1+\\exp{n^{\\prime} } }} )} + \\varepsilon, $$
                                #      \n where \\(K = \\exp(K^{\\prime})\\), \\(n= 1+ \\exp{n^{\\prime}}\\) and \\( \\varepsilon \\sim N(0, \\sigma^{2}) \\).
                                #      \n"),
                                # br(),
                                # br(),
                                # h4("3-5 MLE summary of Transformed Sips"),
                                # tableOutput("mle_table_sips_trans"),
                                # br(),
                                # h4("3-6. Inverse of numerical Hessian of Transformed Sips"),
                                # tableOutput("hessian_num_original_inverse_sips_trans"),
                                )  # end of mainpanel
                     ) # end of sidebar layout
        ), ## end of tab 2

   

        tabPanel("Data",
                 style='padding-left:50px; padding-right:50px; padding-top:40px; padding-bottom:5px',
                 
                 dataTableOutput("table")
        ), ## end of tab 3  
  
       tabPanel("Method",
                
              dashboardPage(
                        dashboardHeader(title = NULL, titleWidth = 0, disable =T
                                        ),
                        sidebar,
                        body
                      ),
              # style='padding-left:0px; padding-right:10px; padding-top:0px; padding-bottom:10px;
              # margin-top:-50px; margin-bottom:0px, margin-left:0px',
              tags$style(type = "text/css", ".container-fluid {padding-left:0px;
                    padding-right:0px;}"),
              tags$style(type = "text/css", ".navbar {margin-bottom: .0px;padding-left:20px}"),
      #        tags$style(type = "text/css", ".container-fluid .navbar-header
#.navbar-brand {margin-left: 0px;}")
       ), ## end of tab 4  

      navbarMenu("More",
               tabPanel("Resources",
                        br(),
                        br(),
                        fluidRow(
                          column("",width=1),
                          column(
                          box(title = "Github", HTML("R Shiny code is available on github. <br/>  <br/>"), width=15,
                              column(12,align="center",
                                              actionButton(inputId='ab1', label=HTML("<img src='github_logo.png' height='40'/> &nbsp; Go to GitHub"), 
                                                         #         icon = icon("github "), 
                                                           style='padding:10px; font-size:160%',
                                                           
                                                                  onclick ="window.open('https://github.com/jgong9/toxicology_app', '_blank')")
                              )
                        , status = "primary",solidHeader = T, collapsible = TRUE,
                        ), width=5)
                        ,
                        
                          
                        column(
                              box(title = "Personal Shiny Server", HTML("This app is accessible through Shiny Server on Raspberry Pi. <br/>  <br/>"), 
                                  width=15,
                                  column(12, align="center",
                                  actionButton(inputId='rasp_button', label=HTML("<img src='raspberry_pi_logo.png' height='40'/> &nbsp; Visit Shiny Server"), 
                                            #   icon = icon("linux"), 
                                               style='padding:10px; font-size:160%',
                                               onclick ="window.open('http://174.99.0.141:9000/toxicology_app', '_blank')")
                                  
                                 ) , status = "danger",solidHeader = T, collapsible = TRUE,
                              ), width=5)
                  ),  
                  fluidRow(
                    column("",width=1),
                    column(
                      box(title = "Contact", HTML("Contact me for any questions or inquiries. <br/>  <br/>"), 
                          width=15,
                          column(12, align="center",
                                 a(actionButton(inputId = "email2", label = HTML("&nbsp; Email Joonho Gong"), 
                                                icon = icon("envelope", lib = "font-awesome"),
                                                style='padding:10px; font-size:160%')
                                 ,
                                 href="mailto:jgong9@ncsu.edu"
                                 )
                          ) , status = "warning",solidHeader = T, collapsible = TRUE,
                      ), width=5)
                  ,  
                  column(
                    box(title = "R Shiny", HTML("This app is available with Shinyapps.io <br/>  <br/>"), width=15,
                        column(12,align="center",
                               actionButton(inputId='ab1', label=HTML("<img src='rstudio_logo3.png' height='40'/> &nbsp; Visit Shinyapps.io"), 
                                            #         icon = icon("github "), 
                                            style='padding:10px; font-size:160%',
                                            
                                            onclick ="window.open('https://jgong9.shinyapps.io/chemistry_app/', '_blank')")
                        )
                        , status = "info",solidHeader = T, collapsible = TRUE,
                    ), width=5)
                  ),
                  style='padding-left:20px; padding-right:10px; padding-top:0px; padding-bottom:5px'
                ),
                tabPanel("About",
                         fluidPage(
                           includeMarkdown("README.md")
                         ),
                         style='padding-left:30px; padding-right:30px; padding-top:10px; padding-bottom:5px',

                )
         ) ## end of tab 4

    ) # end of navbarPage
) # end of tagList





server <- function(input, output) {
    

    inputVar <- reactive({
        if(input$main =="Glyphosate") {
            k_val <- input$k
        } else if(input$main =="Paraquat"){
            k_val <- input$k2
        }
        value <- c(k_val, input$main, 
                   input$conf_level, input$Rfunction, input$method,
                   input$model)
    })
    
    output$otherString <- renderText({
        value <- inputVar()

        paste0(value[2], ", ", value[1])
    })
 
  
    output$trajectoryPlot <- renderPlot({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      ### 1. Original fit
      fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
                           data=data_kth,
                           start=list(Qmax.mol.kg = initial_values[[1]],
                                      Kd.l.mol = initial_values[[2]]))
      
      est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
      max_y <- max( est_curve_origin )
      summary_1 <-summary(fit_nonBoot)
      
      mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      
      ### 3. MLE
      # 3-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                  sqrt(sigma2_initial))

  
      # 3-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      # 3-3 get estimates
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]

      # 3-4 prepare for plotting
      
      Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * data.matrix(data_kth[,1]) ) / (1 + exp(coef_mle[2]) * data.matrix(data_kth[,1]))
      res_mle <- data_kth[,2] - Y_original_pred_mle
      mse_3 <- (sum(res_mle^(2)) ) / (n_data_kth  )
      # compute mse
      
      
      fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_matrix[kth, ]) / (1 + exp(coef_mle[2]) * x_matrix[kth, ])
      # fitted values based on mle for plotting 
      
      
      # 3-5 MLE Confidence band 
      if( input$conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est <- solve(fit_mle$hessian)
      c_pred_vec <- x_matrix[kth, ]
      
      mle_band_upper <- rep(NA, length(c_pred_vec))
      mle_band_lower <- rep(NA, length(c_pred_vec))
      
      for(i in 1:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                sigma = sigma_est_mle, cov_mat=Cov_est)
        
        term <- conf_val * se_est
        mle_band_upper[i] <- fitted_values_mle[i] + term
        mle_band_lower[i] <- fitted_values_mle[i] - term
      }
      
      
      ### 4. Plot
      
      max_final <- max(c(max_y, Y_original_pred, unlist(data_kth[-1,2]), mle_band_upper ))
      
      # main plot
      plot(x_matrix[kth, ], est_curve_origin, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
           ylim = c(0, max_final ), lwd = 2, col= 'blue',
           main=paste("NonTrans-NonLinear-OLS vs Trans-NonLinear-MLE; ", value[2], ", Data ID: ", value[1]))
      lines(c(x_com_vec), c(Y_com_pred), col="green", lty=2, lwd=2)
      lines(x_matrix[kth, ], fitted_values_mle, col="red", lty=3, lwd=2)
      
      # Confidence band
      lines(c_pred_vec, mle_band_upper, col="black", lty=3, lwd=2)
      lines(c_pred_vec, mle_band_lower, col="black", lty=3, lwd=2)
      
      # observed data points
      points(data_set[,(2* kth -1):(2*kth)])
      
      # legend
      legend("bottomright",legend=c(
        paste("NonTrans+NonLinear+OLS : Qm=",
              signif(coefficients(fit_nonBoot)["Qmax.mol.kg"],6),
              ", Kd=",round(coefficients(fit_nonBoot)["Kd.l.mol"]),
              ", MSE=",format(mse_1, scientific = T),sep=""),
        paste("Trans+NonLinear+MLE : Qm=",
              signif(Q_est_mle,6),
              ", Kd=",round(K_est_mle),
              ", MSE=",format(mse_3, scientific = T),sep=""),
        paste("Commercial Product : Qm=",
              signif(Q_com,6),
              ", Kd=",round(K_com),
              ", MSE=",format(mse_0,scientific = T),sep=""),
        paste( input$conf_level, " Confidence band of MLE"),
        paste("Observed data point")
      ),
      bty='n',
      col=c("blue","red","green", "black", "black"),lwd=c(2, 2, 2, 2, NA),lty=c(1,2,3, 3,NA),
      density = c(0, 0, 0, 0, NA), fill = c("blue", "red","green", "black", "white"),
      border=c(NA, NA, NA, NA, NA),
      pch = c(NA, NA, NA, NA, 1),
      cex=1
      )
    })
    
    
    
    output$mle_table <- renderTable({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      ### 1. Original fit
      fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
                           data=data_kth,
                           start=list(Qmax.mol.kg = initial_values[[1]],
                                      Kd.l.mol = initial_values[[2]]))
      
      est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
      max_y <- max( est_curve_origin )
      summary_1 <-summary(fit_nonBoot)
      
      mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute rse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      ### 3. MLE
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
  
      
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       #    ndigit = 20,
                       #    method = "CG",
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      
      
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]

      
      Cov_est <- solve(fit_mle$hessian)
      std_error1 <- sqrt(Cov_est[1,1]  )
      std_error2 <- sqrt(Cov_est[2,2]  )
      std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
      std_error3 <- sqrt(Cov_est[3,3]  )
      
      if( input$conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if (input$conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
      SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
      upper <- c(coef_mle[1] + conf_val * std_error1,
                 coef_mle[2] + conf_val * std_error2, 
                 exp(coef_mle[2] + conf_val * std_error2),
                 # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                 coef_mle[3] + conf_val * std_error3 )
      lower <- c(coef_mle[1] - conf_val * std_error1,
                 coef_mle[2] - conf_val * std_error2, 
                 exp(coef_mle[2] - conf_val * std_error2),
                 # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                 coef_mle[3] - conf_val * std_error3 )
      
      

      table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T), Lower=format(lower, scientific = T))
      

    },
    bordered=T)
    
    
    output$aic_lan <- renderText({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      ### 1. Original fit
      fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
                           data=data_kth,
                           start=list(Qmax.mol.kg = initial_values[[1]],
                                      Kd.l.mol = initial_values[[2]]))
      
      est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
      max_y <- max( est_curve_origin )
      summary_1 <-summary(fit_nonBoot)
      
      mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute rse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      ### 3. MLE
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       #    ndigit = 20,
                       #    method = "CG",
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      
      
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]
      

      
      AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2]))
      
      paste("AIC of Langmuir model is ", AIC_lan)
          
      })
    
    
    
    # output$hessian_num_original <- renderTable({
    #   value <- inputVar()
    #   
    #   kth <- as.numeric(value[1])
    #   
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    #   
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    #   
    #   
    #   ### 0. Commerical product
    #   Q_com <- initial_values[[1]]
    #   K_com <- initial_values[[2]]
    #   
    #   x_com_vec <- x_matrix[kth,]
    #   Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
    #   
    #   mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
    #   
    #   ### 1. Original fit
    #   fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
    #                        data=data_kth,
    #                        start=list(Qmax.mol.kg = initial_values[[1]],
    #                                   Kd.l.mol = initial_values[[2]]))
    #   
    #   est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
    #   max_y <- max( est_curve_origin )
    #   summary_1 <-summary(fit_nonBoot)
    #   
    #   mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
    #   
    #   
    #   ### 2. Linear model
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    #   
    #   data_trans <- data.frame(q.mol.kg = Y_trans, 
    #                            Cw.mol.l = X_trans)
    #   
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
    #   
    #   x_vec <- 1 / x_matrix[kth,-1]
    #   Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
    #   Y_original_pred <- 1/ Y_trans_pred
    #   
    #   # Compute rse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    #   
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    #   
    #   ### 3. MLE
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
    #                           sqrt(sigma2_initial))
    # 
    #   
    #   if(value[4] == "optim"){
    #     fit_mle <- optim(fn = negLogLikelihood,
    #                      par = initial_value_mle,
    #                      hessian = T,
    #                      method = value[5],
    #                      data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle <- fit_mle$par
    #     
    #   } else if(value[4] == "nlm"){
    #     fit_mle <- nlm(f = negLogLikelihood,
    #                    p = initial_value_mle,
    #                    hessian = T,
    #                    data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     
    #     coef_mle <- fit_mle$estimate
    #     
    #   }
    #   
    #   
    # 
    #   Q_est_mle <- coef_mle[1]
    #   K_est_mle <- exp( coef_mle[2] )
    #   sigma_est_mle <- coef_mle[3]
    # 
    #   
    #   fit_mle$hessian
    #   
    #   
    #   
    # },
    # colnames = F,
    # digits=-2,
    # bordered =T
    # )
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # output$hessian_num_original_inverse <- renderTable({
    #   value <- inputVar()
    #   
    #   kth <- as.numeric(value[1])
    #   
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    #   
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    #   
    #   
    #   ### 0. Commerical product
    #   Q_com <- initial_values[[1]]
    #   K_com <- initial_values[[2]]
    #   
    #   x_com_vec <- x_matrix[kth,]
    #   Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
    #   
    #   mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
    #   
    #   ### 1. Original fit
    #   fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
    #                        data=data_kth,
    #                        start=list(Qmax.mol.kg = initial_values[[1]],
    #                                   Kd.l.mol = initial_values[[2]]))
    #   
    #   est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
    #   max_y <- max( est_curve_origin )
    #   summary_1 <-summary(fit_nonBoot)
    #   
    #   mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
    #   
    #   
    #   ### 2. Linear model
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    #   
    #   data_trans <- data.frame(q.mol.kg = Y_trans, 
    #                            Cw.mol.l = X_trans)
    #   
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
    #   
    #   x_vec <- 1 / x_matrix[kth,-1]
    #   Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
    #   Y_original_pred <- 1/ Y_trans_pred
    #   
    #   # Compute rse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    #   
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    #   
    #   ### 3. MLE
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
    #                           sqrt(sigma2_initial))
    # 
    #   
    #   if(value[4] == "optim"){
    #     fit_mle <- optim(fn = negLogLikelihood,
    #                      par = initial_value_mle,
    #                      hessian = T,
    #                      method = value[5],
    #                      data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle <- fit_mle$par
    #     
    #   } else if(value[4] == "nlm"){
    #     fit_mle <- nlm(f = negLogLikelihood,
    #                    p = initial_value_mle,
    #                    hessian = T,
    #                    data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     
    #     coef_mle <- fit_mle$estimate
    #     
    #   }
    #   
    #   
    # 
    #   Q_est_mle <- coef_mle[1]
    #   K_est_mle <- exp( coef_mle[2] )
    #   sigma_est_mle <- coef_mle[3]
    # 
    #   solve(fit_mle$hessian)
    #   
    # },
    # colnames = F,
    # digits=-2,
    # bordered = T
    # )
    
###################################################################
#### Model 2. Freundlich
    
    
    output$plot_fre <- renderPlot({
      value <- inputVar()
      
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      
      ### 3. MLE
      # 3-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      # 3-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle_lan <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle_lan <- fit_mle_lan$par
        
      } else if(value[4] == "nlm"){
        fit_mle_lan <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle_lan <- fit_mle_lan$estimate
        
      }
      
      # 3-3 get estimates
      Q_est_mle_lan <- coef_mle_lan[1]
      K_est_mle_lan <- exp( coef_mle_lan[2] )
      sigma_est_mle_lan <- coef_mle_lan[3]
      
      # 3-4 prepare for plotting
      
      Y_original_pred_mle_lan <- ( coef_mle_lan[1] * exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]) ) / (1 + exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]))
      res_mle_lan <- data_kth[,2] - Y_original_pred_mle_lan
      mse_lan <- (sum(res_mle_lan^(2)) ) / (n_data_kth  )
      # compute mse
      
      
      fitted_values_mle_lan <-  ( coef_mle_lan[1] * exp(coef_mle_lan[2]) *  x_matrix[kth, ]) / (1 + exp(coef_mle_lan[2]) * x_matrix[kth, ])
      # fitted values based on mle for plotting 
      
      
      # 3-5 MLE Confidence band 
      if( input$conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est_lan <- solve(fit_mle_lan$hessian)
      c_pred_vec <- x_matrix[kth, ]
      
      mle_band_upper_lan <- rep(NA, length(c_pred_vec))
      mle_band_lower_lan <- rep(NA, length(c_pred_vec))
      
      for(i in 1:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle(Q = Q_est_mle_lan, K_new=coef_mle_lan[2], C=c_now,
                                sigma = sigma_est_mle_lan, cov_mat=Cov_est_lan
                                )
        
        term <- conf_val * se_est
        mle_band_upper_lan[i] <- fitted_values_mle_lan[i] + term
        mle_band_lower_lan[i] <- fitted_values_mle_lan[i] - term
      }
      
      
      
      
      
      
      
      
      
      
      ### 2. Freundlich model
      ## 2-1 initial value
      Y_trans <-  log(data_kth[-1,2] )
      X_trans <- log(data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      
      
      res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
        # ## 2-2. MLE
        # initial_value_mle <- c( K_F_est_lm, n_est_lm, sqrt(sigma2_initial))
        # 
        # 
        # # 2-2 fit a new model of reparameterization for Kd
        # if(value[4] == "optim"){
        #   fit_mle_fre <- optim(fn = negLogLikelihood_fre,
        #                    par = initial_value_mle,
        #                    hessian = T,
        #                    method = value[5],
        #                    data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        #   )
        #   coef_mle_fre <- fit_mle_fre$par
        #   
        # } else if(value[4] == "nlm"){
        #   fit_mle_fre <- nlm(f = negLogLikelihood_fre,
        #                  p = initial_value_mle,
        #                  hessian = T,
        #                  data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        #   )
        #   
        #   coef_mle_fre <- fit_mle_fre$estimate
        #   
        # }
        # 
        # # 2-3 get estimates
        # K_F_est_mle_fre <- coef_mle_fre[1]
        # n_est_mle_fre <- coef_mle_fre[2]
        # sigma_est_mle_fre <- coef_mle_fre[3]
        # 
        # # 2-4 prepare for plotting
        # Y_original_pred_mle_fre <-  coef_mle_fre[1] * data.matrix(data_kth[,1])^(1/coef_mle_fre[2]) 
        # res_mle_fre <- data_kth[,2] - Y_original_pred_mle_fre
        # mse_fre <- (sum(res_mle_fre^(2)) ) / (n_data_kth  )
        # # compute mse
        # 
        # 
        # fitted_values_mle_fre <-   coef_mle_fre[1] *   x_matrix[kth, ]^(1/coef_mle_fre[2])
        # # fitted values based on mle for plotting 
        # 
        # 
        # # 2-5 MLE Confidence band 
        # if( input$conf_level == "95%"){
        #   conf_val <- qnorm(0.975)
        # } else if( input$conf_level == "90%"){
        #   conf_val <- qnorm(0.95)
        # }  else if (input$conf_level == "80%"){
        #   conf_val<- qnorm(0.90)
        # } else if(input$conf_level == "99%"){
        #   conf_val<- qnorm(0.995)
        # }
        # 
        # Cov_est_fre <- solve(fit_mle_fre$hessian)
        # c_pred_vec <- x_matrix[kth, ]
        # 
        # mle_band_upper_fre <- rep(NA, length(c_pred_vec))
        # mle_band_lower_fre <- rep(NA, length(c_pred_vec))
        # 
        # for(i in 2:length(c_pred_vec)){
        #   c_now <- c_pred_vec[i]
        #   se_est <- se_for_CI_mle_fre( K_F = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, 
        #                                sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre
        #                                )
        #   
        #   term <- conf_val * se_est
        #   mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
        #   mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
        # }
        # 
        # mle_band_upper_fre[1] <- 0
        # mle_band_lower_fre[1] <- 0
        # # since derivative[2] does not exist when C= 0 
      
        ## 2-2. MLE
        initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))


        # 2-2 fit a new model of reparameterization for Kd
        if(value[4] == "optim"){
          fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                           par = initial_value_mle,
                           hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
          )
          coef_mle_fre <- fit_mle_fre$par

        } else if(value[4] == "nlm"){
          fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                         p = initial_value_mle,
                         hessian = T,
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
          )

          coef_mle_fre <- fit_mle_fre$estimate

        }

        # 2-3 get estimates
        K_F_prime_est_mle_fre <- coef_mle_fre[1]
        K_F_est_mle_fre <- exp(coef_mle_fre[1])

        n_est_mle_fre <- coef_mle_fre[2]
        sigma_est_mle_fre <- coef_mle_fre[3]

        # 2-4 prepare for plotting
        Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * data.matrix(data_kth[,1])^(1/coef_mle_fre[2])
        res_mle_fre <- data_kth[,2] - Y_original_pred_mle_fre
        mse_fre <- (sum(res_mle_fre^(2)) ) / (n_data_kth  )
        # compute mse


        fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_matrix[kth, ]^(1/coef_mle_fre[2])
        # fitted values based on mle for plotting


        # 2-5 MLE Confidence band
        if( input$conf_level == "95%"){
          conf_val <- qnorm(0.975)
        } else if( input$conf_level == "90%"){
          conf_val <- qnorm(0.95)
        }  else if (input$conf_level == "80%"){
          conf_val<- qnorm(0.90)
        } else if(input$conf_level == "99%"){
          conf_val<- qnorm(0.995)
        }

        Cov_est_fre <- solve(fit_mle_fre$hessian)
        c_pred_vec <- x_matrix[kth, ]

        mle_band_upper_fre <- rep(NA, length(c_pred_vec))
        mle_band_lower_fre <- rep(NA, length(c_pred_vec))

        for(i in 2:length(c_pred_vec)){
          c_now <- c_pred_vec[i]
          se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)

          term <- conf_val * se_est
          mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
          mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
        }

        mle_band_upper_fre[1] <- 0
        mle_band_lower_fre[1] <- 0
        # since derivative[2] does not exist when C= 0

      ### 4. Plot
      
      max_final <- max(c(unlist(data_kth[-1,2]), mle_band_upper_fre, mle_band_upper_lan ))
      
      # main plot
      plot(x_matrix[kth, ], fitted_values_mle_fre, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
           ylim = c(0, max_final ), lwd = 2, col= 'blue',
           main=paste( value[2], ", Data ID: ", value[1]))
      
      
      # Confidence band
      lines(c_pred_vec, mle_band_upper_fre, col="blue", lty=3, lwd=2)
      lines(c_pred_vec, mle_band_lower_fre, col="blue", lty=3, lwd=2)
      
      
      lines(c_pred_vec, fitted_values_mle_lan, col="red", lty=1, lwd=2)
      
      lines(c_pred_vec, mle_band_upper_lan, col="red", lty=3, lwd=2)
      lines(c_pred_vec, mle_band_lower_lan, col="red", lty=3, lwd=2)
      lines(c(x_com_vec), c(Y_com_pred), col="green", lty=2, lwd=2)
      
      # observed data points
      points(data_set[,(2* kth -1):(2*kth)])
      
  
      legend("bottomright",legend=c(
        paste("Langmuir : Qm=",
              signif(Q_est_mle_lan,6),
              ", Kd=",round(K_est_mle_lan),
              ", MSE=",format(mse_lan,scientific = T),sep=""),
        paste("Freundlich : K_F=",
              signif(K_F_est_mle_fre),
              ", n=",signif(n_est_mle_fre),
              ", MSE=", format(mse_fre,scientific = T ),sep=""),
        paste( input$conf_level, " Confidence band of Langmuir MLE"),
        paste( input$conf_level, " Confidence band of Freundlich MLE"),
        paste("Commercial Product : Qm=",
              signif(Q_com,6),
              ", Kd=",round(K_com),
              ", MSE=",format(mse_0, scientific = T),sep=""),
        paste("Observed data point")
      ),
      bty='n',
      col=c("red","blue","red", "blue", "green", "black"),lwd=c(2, 2, 2, 2,2, NA),lty=c(1,1,3,3,2,NA),
      density = c(0, 0, 0, 0, 0, NA), fill = c("red", "blue","red", "blue", "green", "white"),
      border=c(NA, NA, NA, NA, NA, NA),
      pch = c(NA, NA, NA, NA, NA, 1),
      cex=1
      )
    })
    
    
    
    output$mle_table_fre <- renderTable({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product

      
      ### 2. Freundlich model
      ## 2-1 initial value
      Y_trans <-  log(data_kth[-1,2] )
      X_trans <- log(data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      
      res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      

      
      # 2-2. MLE
       initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))


       # 2-2 fit a new model of reparameterization for Kd
       if(value[4] == "optim"){
         fit_mle <- optim(fn = negLogLikelihood_fre_trans,
                          par = initial_value_mle,
                          hessian = T,
                          method = value[5],
                          data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
         )
         coef_mle <- fit_mle$par

       } else if(value[4] == "nlm"){
         fit_mle <- nlm(f = negLogLikelihood_fre_trans,
                        p = initial_value_mle,
                        hessian = T,
                        data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
         )

         coef_mle <- fit_mle$estimate

       }



       K_F_prime_est_mle_fre <- coef_mle[1]
       K_F_est_mle_fre <- exp(coef_mle[1])

       n_est_mle_fre <- coef_mle[2]
       sigma_est_mle <- coef_mle[3]


       Cov_est <- solve(fit_mle$hessian)
       std_error1 <- sqrt(Cov_est[1,1])
       std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle[1]))^(2))
       std_error2 <- sqrt(Cov_est[2,2])
       std_error3 <- sqrt(Cov_est[3,3])

       if( input$conf_level == "95%"){
         conf_val <- qnorm(0.975)
       } else if( input$conf_level == "90%"){
         conf_val <- qnorm(0.95)
       }  else if (input$conf_level == "80%"){
         conf_val<- qnorm(0.90)
       } else if (input$conf_level == "99%"){
         conf_val<- qnorm(0.995)
       }

       MLE <- c(coef_mle[1], exp(coef_mle[1]), coef_mle[2], coef_mle[3])
       SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
       upper <- c(coef_mle[1] + conf_val * std_error1, 
                  # exp(coef_mle[1]) + conf_val * std_error1_delta_method,
                  exp(coef_mle[1] + conf_val * std_error1),
                  coef_mle[2] + conf_val * std_error2, 
                  coef_mle[3] + conf_val * std_error3 )
       lower <- c(coef_mle[1] - conf_val * std_error1, 
                  # exp(coef_mle[1]) - conf_val * std_error1_delta_method,
                  exp(coef_mle[1] - conf_val * std_error1),
                  coef_mle[2] - conf_val * std_error2, 
                  coef_mle[3] - conf_val * std_error3 )


      table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                 MLE = format(MLE, scientific = T),
                                 se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T),
                                 Lower=format(lower, scientific = T))

    },
    bordered=T)
    
    
    output$aic_fre <- renderText({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      
      
      ### 2. Freundlich model
      ## 2-1 initial value
      Y_trans <-  log(data_kth[-1,2] )
      X_trans <- log(data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      
      res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
      
      
      # 2-2. MLE
      initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
      
      
      # 2-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood_fre_trans,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood_fre_trans,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      
      K_F_prime_est_mle_fre <- coef_mle[1]
      K_F_est_mle_fre <- exp(coef_mle[1])
      
      n_est_mle_fre <- coef_mle[2]
      sigma_est_mle <- coef_mle[3]
      
      
      ### Langmuir
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      # 3-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      # 3-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle_lan <- optim(fn = negLogLikelihood,
                             par = initial_value_mle,
                             hessian = T,
                             method = value[5],
                             data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle_lan <- fit_mle_lan$par
        
      } else if(value[4] == "nlm"){
        fit_mle_lan <- nlm(f = negLogLikelihood,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle_lan <- fit_mle_lan$estimate
        
      }

      
      
      
      AIC_fre <- 2 * length(coef_mle) + negLogLikelihood_fre_trans(coef_mle, data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2]))
      AIC_lan <- 2 * length(coef_mle_lan) + negLogLikelihood(coef_mle_lan, data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2]))
      if(AIC_lan < AIC_fre){
        AIC_winner <- "Langmuir"
      } else if(AIC_lan > AIC_fre){
        AIC_winner <- "Freundlich"
      } else { AIC_winner <- "even" }
      
      
      paste("AIC of Freundlich model is ", AIC_fre, "\nAIC of Langmuir model is ", AIC_lan, "\nWinner is ", AIC_winner)
      
    })
    
    
    
    # output$hessian_num_original_fre <- renderTable({
    #   value <- inputVar()
    #   
    # 
    #   kth <- as.numeric(value[1])
    #   
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    #   
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    #   
    #   
    #   ### 0. Commerical product
    #   
    #   
    #   ### 2. Freundlich model
    #   ## 2-1 initial value
    #   Y_trans <-  log(data_kth[-1,2] )
    #   X_trans <- log(data_kth[-1,1] )
    #   
    #   data_trans <- data.frame(q.mol.kg = Y_trans, 
    #                            Cw.mol.l = X_trans)
    #   
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   K_F_est_lm <- exp( coef_lm[1] )
    #   n_est_lm <- ( 1 / coef_lm[2] ) 
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    #   
    #   
    #   res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   
    #   ## 2-2. MLE
    #   initial_value_mle <- c( K_F_est_lm, n_est_lm, sqrt(sigma2_initial))
    #   
    #   
    #   # 2-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle <- optim(fn = negLogLikelihood_fre,
    #                      par = initial_value_mle,
    #                      hessian = T,
    #                      method = value[5],
    #                      data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle <- fit_mle$par
    #     
    #   } else if(value[4] == "nlm"){
    #     fit_mle <- nlm(f = negLogLikelihood_fre,
    #                    p = initial_value_mle,
    #                    hessian = T,
    #                    data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     
    #     coef_mle <- fit_mle$estimate
    #     
    #   }
    #   
    #   
    # 
    #   
    #   fit_mle$hessian
    #   
    #   
    #   
    # },
    # colnames = F,
    # digits=-2,
    # bordered =T
    # )
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # output$hessian_num_original_inverse_fre <- renderTable({
    #   value <- inputVar()
    #   
    #   kth <- as.numeric(value[1])
    #   
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    #   
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    #   
    #   
    #   ### 0. Commerical product
    #   
    #   
    #   ### 2. Freundlich model
    #   ## 2-1 initial value
    #   Y_trans <-  log(data_kth[-1,2] )
    #   X_trans <- log(data_kth[-1,1] )
    #   
    #   data_trans <- data.frame(q.mol.kg = Y_trans, 
    #                            Cw.mol.l = X_trans)
    #   
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   K_F_est_lm <- exp( coef_lm[1] )
    #   n_est_lm <- ( 1 / coef_lm[2] ) 
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    #   
    #   
    #   res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   
    #   ## 2-2. MLE
    #   initial_value_mle <- c( K_F_est_lm, n_est_lm, sqrt(sigma2_initial))
    #   
    #   
    #   # 2-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle <- optim(fn = negLogLikelihood_fre,
    #                      par = initial_value_mle,
    #                      hessian = T,
    #                      method = value[5],
    #                      data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle <- fit_mle$par
    #     
    #   } else if(value[4] == "nlm"){
    #     fit_mle <- nlm(f = negLogLikelihood_fre,
    #                    p = initial_value_mle,
    #                    hessian = T,
    #                    data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     
    #     coef_mle <- fit_mle$estimate
    #     
    #   }
    # 
    #   
    #   solve(fit_mle$hessian)
    #   
    # },
    # colnames = F,
    # digits=-2,
    # bordered = T
    # )
    
    
    
    
    # ###################################################################
    # #### Model 3. Sips
    # 
    # 
    # 
    # ### Trans
    # output$plot_sips_trans <- renderPlot({
    #   value <- inputVar()
    # 
    # 
    #   kth <- as.numeric(value[1])
    # 
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    # 
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    # 
    # 
    #   ### 0. Commerical product
    #   Q_com <- initial_values[[1]]
    #   K_com <- initial_values[[2]]
    # 
    #   x_com_vec <- x_matrix[kth,]
    #   Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
    # 
    #   mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
    # 
    # 
    #   ### 2. Linear model
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   x_vec <- 1 / x_matrix[kth,-1]
    #   Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
    #   Y_original_pred <- 1/ Y_trans_pred
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 3. MLE
    #   # 3-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]),
    #                           sqrt(sigma2_initial))
    # 
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_lan <- optim(fn = negLogLikelihood,
    #                          par = initial_value_mle,
    #                          hessian = T,
    #                          method = value[5],
    #                          data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_lan <- fit_mle_lan$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_lan <- nlm(f = negLogLikelihood,
    #                        p = initial_value_mle,
    #                        hessian = T,
    #                        data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_lan <- fit_mle_lan$estimate
    # 
    #   }
    # 
    #   # 3-3 get estimates
    #   Q_est_mle_lan <- coef_mle_lan[1]
    #   K_est_mle_lan <- exp( coef_mle_lan[2] )
    #   sigma_est_mle_lan <- coef_mle_lan[3]
    # 
    #   # 3-4 prepare for plotting
    # 
    #   Y_original_pred_mle_lan <- ( coef_mle_lan[1] * exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]) ) / (1 + exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]))
    #   res_mle_lan <- data_kth[,2] - Y_original_pred_mle_lan
    #   mse_lan <- (sum(res_mle_lan^(2)) ) / (n_data_kth  )
    #   # compute mse
    # 
    # 
    #   fitted_values_mle_lan <-  ( coef_mle_lan[1] * exp(coef_mle_lan[2]) *  x_matrix[kth, ]) / (1 + exp(coef_mle_lan[2]) * x_matrix[kth, ])
    #   # fitted values based on mle for plotting
    # 
    # 
    #   # 3-5 MLE Confidence band
    #   if( input$conf_level == "95%"){
    #     conf_val <- qnorm(0.975)
    #   } else if( input$conf_level == "90%"){
    #     conf_val <- qnorm(0.95)
    #   }  else if (input$conf_level == "80%"){
    #     conf_val<- qnorm(0.90)
    #   } else if(input$conf_level == "99%"){
    #     conf_val<- qnorm(0.995)
    #   }
    # 
    #   Cov_est_lan <- solve(fit_mle_lan$hessian)
    #   c_pred_vec <- x_matrix[kth, ]
    # 
    #   mle_band_upper_lan <- rep(NA, length(c_pred_vec))
    #   mle_band_lower_lan <- rep(NA, length(c_pred_vec))
    # 
    #   for(i in 1:length(c_pred_vec)){
    #     c_now <- c_pred_vec[i]
    #     se_est <- se_for_CI_mle(Q = Q_est_mle_lan, K_new=coef_mle_lan[2], C=c_now,
    #                             sigma = sigma_est_mle_lan, cov_mat=Cov_est_lan
    #     )
    # 
    #     term <- conf_val * se_est
    #     mle_band_upper_lan[i] <- fitted_values_mle_lan[i] + term
    #     mle_band_lower_lan[i] <- fitted_values_mle_lan[i] - term
    #   }
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    #   ### Sips model
    #   ## 1. initial value
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 2. MLE
    #   # 2-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]),
    #                           1,
    #                           sqrt(sigma2_initial))
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_sips <- optim(fn = negLogLikelihood_sips_trans,
    #                           par = initial_value_mle,
    #                           hessian = T,
    #                           method = value[5],
    #                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_sips <- fit_mle_sips$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_sips <- nlm(f = negLogLikelihood_sips_trans,
    #                         p = initial_value_mle,
    #                         hessian = T,
    #                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_sips <- fit_mle_sips$estimate
    # 
    #   }
    # 
    #   # 3-3 get estimates
    #   Q_est_mle_sips <- coef_mle_sips[1]
    #   K_est_mle_sips <- exp( coef_mle_sips[2] )
    #   n_est_mle_sips <- 1 + exp( coef_mle_sips[3] )
    #   # n_est_mle_sips <- exp( coef_mle_sips[3] )
    #   ## without 1+
    #   sigma_est_mle_sips <- coef_mle_sips[4]
    # 
    #   # 3-4 prepare for plotting
    # 
    #   Y_original_pred_mle_sips <- ( coef_mle_sips[1] * exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(1/ (1+exp(coef_mle_sips[3]) )  ) ) / (1 + exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(1/ (1+exp(coef_mle_sips[3]) )  ) )
    #   # Y_original_pred_mle_sips <- ( coef_mle_sips[1] * exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(1/ (exp(coef_mle_sips[3]) )  ) ) / (1 + exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(1/ (exp(coef_mle_sips[3]) )  ) )
    #   ## without 1+
    #   res_mle_sips <- data_kth[,2] - Y_original_pred_mle_sips
    #   mse_sips <- (sum(res_mle_sips^(2)) ) / (n_data_kth  )
    #   # compute mse
    # 
    # 
    #   fitted_values_mle_sips <-  ( coef_mle_sips[1] * exp(coef_mle_sips[2]) *  x_matrix[kth, ]^(1/ (1+exp(coef_mle_sips[3]) )  ) ) / (1 + exp(coef_mle_sips[2]) * x_matrix[kth, ]^(1/ (1+exp(coef_mle_sips[3]) )  )  )
    #   # fitted_values_mle_sips <-  ( coef_mle_sips[1] * exp(coef_mle_sips[2]) *  x_matrix[kth, ]^(1/ (exp(coef_mle_sips[3]) )  ) ) / (1 + exp(coef_mle_sips[2]) * x_matrix[kth, ]^(1/ (exp(coef_mle_sips[3]) )  )  )
    #   ## without 1+
    #   # fitted values based on mle for plotting
    # 
    # 
    #   # 3-5 MLE Confidence band
    #   if( input$conf_level == "95%"){
    #     conf_val <- qnorm(0.975)
    #   } else if( input$conf_level == "90%"){
    #     conf_val <- qnorm(0.95)
    #   }  else if (input$conf_level == "80%"){
    #     conf_val<- qnorm(0.90)
    #   } else if(input$conf_level == "99%"){
    #     conf_val<- qnorm(0.995)
    #   }
    # 
    #   Cov_est_sips <- solve(fit_mle_sips$hessian)
    #   c_pred_vec <- x_matrix[kth, ]
    # 
    #   mle_band_upper_sips <- rep(NA, length(c_pred_vec))
    #   mle_band_lower_sips <- rep(NA, length(c_pred_vec))
    # 
    #   for(i in 1:length(c_pred_vec)){
    #     c_now <- c_pred_vec[i]
    #     se_est <- se_for_CI_mle_sips_trans(Q = Q_est_mle_sips, K_new=coef_mle_sips[2],
    #                                  n_new=coef_mle_sips[3], C=c_now,
    #                                  sigma = sigma_est_mle_sips, cov_mat=Cov_est_sips
    #     )
    # 
    #     term <- conf_val * se_est
    #     mle_band_upper_sips[i] <- fitted_values_mle_sips[i] + term
    #     mle_band_lower_sips[i] <- fitted_values_mle_sips[i] - term
    #   }
    # 
    # 
    # 
    # 
    # 
    #   ### 4. Plot
    #   max_val <- max(mle_band_upper_sips, na.rm = T)
    #   max_final <- max(c(unlist(data_kth[-1,2]), max_val, mle_band_upper_lan ))
    # 
    #   # main plot
    #   plot(x_matrix[kth, ], fitted_values_mle_sips, type='l', xlab = 'Cw.mol.l', ylab = 'q.mol.kg',
    #        ylim = c(0, max_final ), lwd = 2, col= 'blue',
    #        main=paste( value[2], ", Data ID: ", value[1]))
    # 
    # 
    #   # Confidence band
    #   lines(c_pred_vec, mle_band_upper_sips, col="blue", lty=3, lwd=2)
    #   lines(c_pred_vec, mle_band_lower_sips, col="blue", lty=3, lwd=2)
    # 
    # 
    #   lines(c_pred_vec, fitted_values_mle_lan, col="red", lty=1, lwd=2)
    # 
    #   lines(c_pred_vec, mle_band_upper_lan, col="red", lty=3, lwd=2)
    #   lines(c_pred_vec, mle_band_lower_lan, col="red", lty=3, lwd=2)
    #   lines(c(x_com_vec), c(Y_com_pred), col="green", lty=2, lwd=2)
    # 
    #   # observed data points
    #   points(data_set[,(2* kth -1):(2*kth)])
    # 
    # 
    #   legend("bottomright",legend=c(
    #     paste("Langmuir : Qm=",
    #           signif(Q_est_mle_lan,6),
    #           ", Kd=",round(K_est_mle_lan),
    #           ", MSE=",format(mse_lan,scientific = T),sep=""),
    #     paste("Sips : Qn=",
    #           signif(Q_est_mle_sips,6),
    #           ", Kd=",round(K_est_mle_sips),
    #           ", n=", round(n_est_mle_sips),
    #           ", MSE=",format(mse_sips,scientific = T),sep=""),
    #     paste( input$conf_level, " Confidence band of Langmuir MLE"),
    #     paste( input$conf_level, " Confidence band of Sips MLE"),
    #     paste("Commercial Product : Qm=",
    #           signif(Q_com,6),
    #           ", Kd=",round(K_com),
    #           ", MSE=",format(mse_0, scientific = T),sep=""),
    #     paste("observed point")
    #   ),
    #   bty='n',
    #   col=c("red","blue","red", "blue", "green", "black"),lwd=c(2, 2, 2, 2,2, NA),lty=c(1,1,3,3,2,NA),
    #   density = c(0, 0, 0, 0, 0, NA), fill = c("red", "blue","red", "blue", "green", "white"),
    #   border=c(NA, NA, NA, NA, NA, NA),
    #   pch = c(NA, NA, NA, NA, NA, 1),
    #   cex=1
    #   )
    # })
    # 
    # 
    # 
    # output$mle_table_sips_trans <- renderTable({
    #   value <- inputVar()
    # 
    #   kth <- as.numeric(value[1])
    # 
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    # 
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    # 
    # 
    #   ### 0. Commerical product
    # 
    # 
    #   ### 2. Freundlich model
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   # ### 2. MLE
    #   # # 2-1 initial values from the previous linear model
    #   # sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   # initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2],),
    #   #                         1,
    #   #                         sqrt(sigma2_initial))
    #   #
    #   # # 3-2 fit a new model of reparameterization for Kd
    #   # if(value[4] == "optim"){
    #   #   fit_mle_sips <- optim(fn = negLogLikelihood_sips,
    #   #                         par = initial_value_mle,
    #   #                         hessian = T,
    #   #                         method = value[5],
    #   #                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #   #   )
    #   #   coef_mle_sips <- fit_mle_sips$par
    #   #
    #   # } else if(value[4] == "nlm"){
    #   #   fit_mle_sips <- nlm(f = negLogLikelihood_sips,
    #   #                       p = initial_value_mle,
    #   #                       hessian = T,
    #   #                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #   #   )
    #   #
    #   #   coef_mle_sips <- fit_mle_sips$estimate
    #   #
    #   # }
    #   #
    #   # # 3-3 get estimates
    #   # Q_est_mle_sips <- coef_mle_sips[1]
    #   # K_est_mle_sips <- exp( coef_mle_sips[2] )
    #   # n_est_mle_sips <- coef_mle_sips[3]
    #   # sigma_est_mle_sips <- coef_mle_sips[4]
    # 
    #   # 2-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]),
    #                           1,
    #                           sqrt(sigma2_initial))
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_sips <- optim(fn = negLogLikelihood_sips_trans,
    #                           par = initial_value_mle,
    #                           hessian = T,
    #                           method = value[5],
    #                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_sips <- fit_mle_sips$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_sips <- nlm(f = negLogLikelihood_sips_trans,
    #                         p = initial_value_mle,
    #                         hessian = T,
    #                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_sips <- fit_mle_sips$estimate
    # 
    #   }
    # 
    # 
    #   Cov_est <- solve(fit_mle_sips$hessian)
    #   std_error1 <- sqrt(Cov_est[1,1]  )
    #   std_error2 <- sqrt(Cov_est[2,2]  )
    #   std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle_sips[2]))^(2) )
    #   std_error3 <- sqrt(Cov_est[3,3])
    #   std_error3_delta_method <- sqrt(Cov_est[3,3] *(exp(coef_mle_sips[3]))^(2)  )
    # 
    #   std_error4 <- sqrt(Cov_est[4,4]  )
    # 
    # 
    #   if( input$conf_level == "95%"){
    #     conf_val <- qnorm(0.975)
    #   } else if( input$conf_level == "90%"){
    #     conf_val <- qnorm(0.95)
    #   }  else if (input$conf_level == "80%"){
    #     conf_val<- qnorm(0.90)
    #   } else if (input$conf_level == "99%"){
    #     conf_val<- qnorm(0.995)
    #   }
    # 
    #   MLE <- c(coef_mle_sips[1], coef_mle_sips[2], exp(coef_mle_sips[2]),
    #            coef_mle_sips[3],
    #            1+ exp(coef_mle_sips[3]),
    #            ##
    #            # exp(coef_mle_sips[3]),
    #            ## without 1+
    #            coef_mle_sips[4])
    #   SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3, std_error3_delta_method,std_error4)
    #   upper <- c(coef_mle_sips[1] + conf_val * std_error1,
    #              coef_mle_sips[2] + conf_val * std_error2,
    #              exp(coef_mle_sips[2] + conf_val * std_error2),
    #              # exp(coef_mle_sips[2]) + conf_val * std_error2_delta_method,
    #              coef_mle_sips[3] + conf_val * std_error3,
    #               1 + exp(coef_mle_sips[3] + conf_val * std_error3),
    #              # exp(coef_mle_sips[3] + conf_val * std_error3),
    #              ## without 1+
    #              coef_mle_sips[4] + conf_val * std_error4
    #   )
    #   lower <- c(coef_mle_sips[1] - conf_val * std_error1,
    #              coef_mle_sips[2] - conf_val * std_error2,
    #              exp(coef_mle_sips[2] - conf_val * std_error2),
    #              # exp(coef_mle_sips[2]) - conf_val * std_error2_delta_method,
    #              coef_mle_sips[3] - conf_val * std_error3,
    #              1+exp(coef_mle_sips[3] - conf_val * std_error3),
    #              # exp(coef_mle_sips[3] - conf_val * std_error3),
    #              ## without 1+
    #              coef_mle_sips[4] - conf_val * std_error4
    #   )
    # 
    # 
    # 
    #   table_result <- data.frame(Par = c("Qmax", "Kprime", "K", "nprime", "n","sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
    #                              Upper=format(upper, scientific = T), Lower=format(lower, scientific = T))
    # 
    # 
    # },
    # bordered=T)
    # 
    # 
    # 
    # 
    # 
    # output$hessian_num_original_inverse_sips_trans <- renderTable({
    #   value <- inputVar()
    # 
    #   kth <- as.numeric(value[1])
    # 
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    # 
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    # 
    # 
    #   ### 0. Commerical product
    # 
    # 
    #   ### 2. Freundlich model
    #   ## 2-1 initial value
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 2. MLE
    #   # 2-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]),
    #                           1,
    #                           sqrt(sigma2_initial))
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_sips <- optim(fn = negLogLikelihood_sips_trans,
    #                           par = initial_value_mle,
    #                           hessian = T,
    #                           method = value[5],
    #                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_sips <- fit_mle_sips$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_sips <- nlm(f = negLogLikelihood_sips_trans,
    #                         p = initial_value_mle,
    #                         hessian = T,
    #                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_sips <- fit_mle_sips$estimate
    # 
    #   }
    # 
    # 
    # 
    #   solve(fit_mle_sips$hessian)
    # 
    # },
    # colnames = F,
    # digits=-2,
    # bordered = T
    # )
    # 
    # 
    # #### No transformation on n
    # output$plot_sips <- renderPlot({
    #   value <- inputVar()
    # 
    # 
    #   kth <- as.numeric(value[1])
    # 
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    # 
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    # 
    # 
    #   ### 0. Commerical product
    #   Q_com <- initial_values[[1]]
    #   K_com <- initial_values[[2]]
    # 
    #   x_com_vec <- x_matrix[kth,]
    #   Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
    # 
    #   mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
    # 
    # 
    #   ### 2. Linear model
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   x_vec <- 1 / x_matrix[kth,-1]
    #   Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
    #   Y_original_pred <- 1/ Y_trans_pred
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 3. MLE
    #   # 3-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]),
    #                           sqrt(sigma2_initial))
    # 
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_lan <- optim(fn = negLogLikelihood,
    #                          par = initial_value_mle,
    #                          hessian = T,
    #                          method = value[5],
    #                          data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_lan <- fit_mle_lan$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_lan <- nlm(f = negLogLikelihood,
    #                        p = initial_value_mle,
    #                        hessian = T,
    #                        data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_lan <- fit_mle_lan$estimate
    # 
    #   }
    # 
    #   # 3-3 get estimates
    #   Q_est_mle_lan <- coef_mle_lan[1]
    #   K_est_mle_lan <- exp( coef_mle_lan[2] )
    #   sigma_est_mle_lan <- coef_mle_lan[3]
    # 
    #   # 3-4 prepare for plotting
    # 
    #   Y_original_pred_mle_lan <- ( coef_mle_lan[1] * exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]) ) / (1 + exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]))
    #   res_mle_lan <- data_kth[,2] - Y_original_pred_mle_lan
    #   mse_lan <- (sum(res_mle_lan^(2)) ) / (n_data_kth  )
    #   # compute mse
    # 
    # 
    #   fitted_values_mle_lan <-  ( coef_mle_lan[1] * exp(coef_mle_lan[2]) *  x_matrix[kth, ]) / (1 + exp(coef_mle_lan[2]) * x_matrix[kth, ])
    #   # fitted values based on mle for plotting
    # 
    # 
    #   # 3-5 MLE Confidence band
    #   if( input$conf_level == "95%"){
    #     conf_val <- qnorm(0.975)
    #   } else if( input$conf_level == "90%"){
    #     conf_val <- qnorm(0.95)
    #   }  else if (input$conf_level == "80%"){
    #     conf_val<- qnorm(0.90)
    #   } else if(input$conf_level == "99%"){
    #     conf_val<- qnorm(0.995)
    #   }
    # 
    #   Cov_est_lan <- solve(fit_mle_lan$hessian)
    #   c_pred_vec <- x_matrix[kth, ]
    # 
    #   mle_band_upper_lan <- rep(NA, length(c_pred_vec))
    #   mle_band_lower_lan <- rep(NA, length(c_pred_vec))
    # 
    #   for(i in 1:length(c_pred_vec)){
    #     c_now <- c_pred_vec[i]
    #     se_est <- se_for_CI_mle(Q = Q_est_mle_lan, K_new=coef_mle_lan[2], C=c_now,
    #                             sigma = sigma_est_mle_lan, cov_mat=Cov_est_lan
    #     )
    # 
    #     term <- conf_val * se_est
    #     mle_band_upper_lan[i] <- fitted_values_mle_lan[i] + term
    #     mle_band_lower_lan[i] <- fitted_values_mle_lan[i] - term
    #   }
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    #   ### Sips model
    #   ## 1. initial value
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 2. MLE
    #   # 2-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], 
    #                           log( coef_lm[1] / coef_lm[2],),
    #                           # ( coef_lm[1] / coef_lm[2]),
    #                           1,
    #                           # 1/2,
    #                           sqrt(sigma2_initial))
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_sips <- optim(fn = negLogLikelihood_sips,
    #                          par = initial_value_mle,
    #                          hessian = T,
    #                          method = value[5],
    #                          data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_sips <- fit_mle_sips$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_sips <- nlm(f = negLogLikelihood_sips,
    #                        p = initial_value_mle,
    #                        hessian = T,
    #                        data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_sips <- fit_mle_sips$estimate
    # 
    #   }
    # 
    #   # 3-3 get estimates
    #   Q_est_mle_sips <- coef_mle_sips[1]
    #   K_est_mle_sips <- exp( coef_mle_sips[2] )
    #   # K_est_mle_sips <- coef_mle_sips[2] 
    #   
    #   n_est_mle_sips <- coef_mle_sips[3]
    #   # n_est_mle_sips <- 1/coef_mle_sips[3]
    #   ## reverse 1/n
    #   sigma_est_mle_sips <- coef_mle_sips[4]
    # 
    #   # 3-4 prepare for plotting
    # 
    #  # Y_original_pred_mle_sips <- ( coef_mle_sips[1] * coef_mle_sips[2] * data.matrix(data_kth[,1])^(1/coef_mle_sips[3]) ) / (1 + coef_mle_sips[2] * data.matrix(data_kth[,1])^(1/coef_mle_sips[3]) )
    #   
    #   Y_original_pred_mle_sips <- ( coef_mle_sips[1] * exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(1/coef_mle_sips[3]) ) / (1 + exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(1/coef_mle_sips[3]) )
    #   # Y_original_pred_mle_sips <- ( coef_mle_sips[1] * exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(coef_mle_sips[3]) ) / (1 + exp(coef_mle_sips[2]) * data.matrix(data_kth[,1])^(coef_mle_sips[3]) )
    #   ## reverse 1/n
    #   res_mle_sips <- data_kth[,2] - Y_original_pred_mle_sips
    #   mse_sips <- (sum(res_mle_sips^(2)) ) / (n_data_kth  )
    #   # compute mse
    # 
    #   # fitted_values_mle_sips <-  ( coef_mle_sips[1] * coef_mle_sips[2] *  x_matrix[kth, ]^(1/coef_mle_sips[3]) ) / (1 + coef_mle_sips[2] * x_matrix[kth, ]^(1/coef_mle_sips[3])  )
    #   
    #  fitted_values_mle_sips <-  ( coef_mle_sips[1] * exp(coef_mle_sips[2]) *  x_matrix[kth, ]^(1/coef_mle_sips[3]) ) / (1 + exp(coef_mle_sips[2]) * x_matrix[kth, ]^(1/coef_mle_sips[3])  )
    #  # fitted_values_mle_sips <-  ( coef_mle_sips[1] * exp(coef_mle_sips[2]) *  x_matrix[kth, ]^(coef_mle_sips[3]) ) / (1 + exp(coef_mle_sips[2]) * x_matrix[kth, ]^(coef_mle_sips[3])  )
    #   ## reverse 1/n
    #   # fitted values based on mle for plotting
    # 
    # 
    #   # 3-5 MLE Confidence band
    #   if( input$conf_level == "95%"){
    #     conf_val <- qnorm(0.975)
    #   } else if( input$conf_level == "90%"){
    #     conf_val <- qnorm(0.95)
    #   }  else if (input$conf_level == "80%"){
    #     conf_val<- qnorm(0.90)
    #   } else if(input$conf_level == "99%"){
    #     conf_val<- qnorm(0.995)
    #   }
    # 
    #   Cov_est_sips <- solve(fit_mle_sips$hessian)
    #   c_pred_vec <- x_matrix[kth, ]
    # 
    #   mle_band_upper_sips <- rep(NA, length(c_pred_vec))
    #   mle_band_lower_sips <- rep(NA, length(c_pred_vec))
    # 
    #   for(i in 1:length(c_pred_vec)){
    #     c_now <- c_pred_vec[i]
    #     se_est <- se_for_CI_mle_sips(Q = Q_est_mle_sips, K_new=coef_mle_sips[2],
    #                                  n=n_est_mle_sips, C=c_now,
    #                             sigma = sigma_est_mle_sips, cov_mat=Cov_est_sips
    #     )
    # 
    #     term <- conf_val * se_est
    #     mle_band_upper_sips[i] <- fitted_values_mle_sips[i] + term
    #     mle_band_lower_sips[i] <- fitted_values_mle_sips[i] - term
    #   }
    # 
    # 
    # 
    # 
    # 
    #   ### 4. Plot
    #   mle_sips_max <- max(mle_band_upper_sips, na.rm = T)
    #   max_final <- max(c(unlist(data_kth[-1,2]),
    #                      mle_sips_max,
    #                      mle_band_upper_lan ))
    # 
    #   # main plot
    #   plot(x_matrix[kth, ], fitted_values_mle_sips, type='l', xlab = 'Cw.mol.l', ylab = 'q.mol.kg',
    #        ylim = c(0, max_final ), lwd = 2, col= 'blue',
    #        main=paste( value[2], ", Data ID: ", value[1]))
    # 
    # 
    #   # Confidence band
    #   lines(c_pred_vec, mle_band_upper_sips, col="blue", lty=3, lwd=2)
    #   lines(c_pred_vec, mle_band_lower_sips, col="blue", lty=3, lwd=2)
    # 
    # 
    #   lines(c_pred_vec, fitted_values_mle_lan, col="red", lty=1, lwd=2)
    # 
    #   lines(c_pred_vec, mle_band_upper_lan, col="red", lty=3, lwd=2)
    #   lines(c_pred_vec, mle_band_lower_lan, col="red", lty=3, lwd=2)
    #   lines(c(x_com_vec), c(Y_com_pred), col="green", lty=2, lwd=2)
    # 
    #   # observed data points
    #   points(data_set[,(2* kth -1):(2*kth)])
    # 
    # 
    #   legend("bottomright",legend=c(
    #     paste("Langmuir : Qm=",
    #           signif(Q_est_mle_lan,6),
    #           ", Kd=",round(K_est_mle_lan),
    #           ", MSE=",format(mse_lan,scientific = T),sep=""),
    #     paste("Sips : Qn=",
    #           signif(Q_est_mle_sips,6),
    #           ", Kd=",signif(K_est_mle_sips,6),
    #           ", n=", signif(n_est_mle_sips,6),
    #           ", MSE=",format(mse_sips,scientific = T),sep=""),
    #     paste( input$conf_level, " Confidence band of Langmuir MLE"),
    #     paste( input$conf_level, " Confidence band of Sips MLE"),
    #     paste("Commercial Product : Qm=",
    #           signif(Q_com,6),
    #           ", Kd=",round(K_com),
    #           ", MSE=",format(mse_0, scientific = T),sep=""),
    #     paste("observed point")
    #   ),
    #   bty='n',
    #   col=c("red","blue","red", "blue", "green", "black"),lwd=c(2, 2, 2, 2,2, NA),lty=c(1,1,3,3,2,NA),
    #   density = c(0, 0, 0, 0, 0, NA), fill = c("red", "blue","red", "blue", "green", "white"),
    #   border=c(NA, NA, NA, NA, NA, NA),
    #   pch = c(NA, NA, NA, NA, NA, 1),
    #   cex=1
    #   )
    # })
    # 
    # 
    # 
    # output$mle_table_sips <- renderTable({
    #   value <- inputVar()
    # 
    #   kth <- as.numeric(value[1])
    # 
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    # 
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    # 
    # 
    #   ### 0. Commerical product
    # 
    # 
    #   ### 2. Freundlich model
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 2. MLE
    #   # 2-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], 
    #                           log( coef_lm[1] / coef_lm[2],),
    #                           #( coef_lm[1] / coef_lm[2] ),
    #                           
    #                           1,
    #                           #1/2,
    #                           sqrt(sigma2_initial))
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_sips <- optim(fn = negLogLikelihood_sips,
    #                           par = initial_value_mle,
    #                           hessian = T,
    #                           method = value[5],
    #                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_sips <- fit_mle_sips$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_sips <- nlm(f = negLogLikelihood_sips,
    #                         p = initial_value_mle,
    #                         hessian = T,
    #                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_sips <- fit_mle_sips$estimate
    # 
    #   }
    # 
    #   # 3-3 get estimates
    #   Q_est_mle_sips <- coef_mle_sips[1]
    #   # K_est_mle_sips <-  coef_mle_sips[2] 
    #   
    #   K_est_mle_sips <- exp( coef_mle_sips[2] )
    #   n_est_mle_sips <- coef_mle_sips[3]
    #   # n_est_mle_sips <- 1/ coef_mle_sips[3]
    #   ## reverse 1/n
    #   sigma_est_mle_sips <- coef_mle_sips[4]
    # 
    # 
    #   Cov_est <- solve(fit_mle_sips$hessian)
    #   std_error1 <- sqrt(Cov_est[1,1]  )
    #   std_error2 <- sqrt(Cov_est[2,2]  )
    #   std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle_sips[2]))^(2) )
    #   std_error3 <- sqrt(Cov_est[3,3])
    #  #  std_error3_reverse <- sqrt(Cov_est[3,3] * (-1 * (coef_mle[3])^(-2))^(2)  )
    #   
    #   std_error4 <- sqrt(Cov_est[4,4]  )
    # 
    # 
    #   if( input$conf_level == "95%"){
    #     conf_val <- qnorm(0.975)
    #   } else if( input$conf_level == "90%"){
    #     conf_val <- qnorm(0.95)
    #   }  else if (input$conf_level == "80%"){
    #     conf_val<- qnorm(0.90)
    #   } else if (input$conf_level == "99%"){
    #     conf_val<- qnorm(0.995)
    #   }
    # 
    #   MLE <- c(coef_mle_sips[1], 
    #            coef_mle_sips[2], 
    #            exp(coef_mle_sips[2]), 
    #            coef_mle_sips[3],
    #            # coef_mle_sips[3], 
    #           #  1/coef_mle_sips[3], 
    #            ## revserse 1/n
    #            coef_mle_sips[4])
    #   SE <- c(std_error1, 
    #           std_error2, 
    #           std_error2_delta_method,
    #           std_error3,
    #           # std_error3_reverse,
    #           std_error4)
    #   upper <- c(coef_mle_sips[1] + conf_val * std_error1,
    #              coef_mle_sips[2] + conf_val * std_error2,
    #              exp(coef_mle_sips[2] + conf_val * std_error2),
    #              # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
    #              coef_mle_sips[3] + conf_val * std_error2,
    #              
    #              # coef_mle_sips[3] + conf_val * std_error3,
    #             #  1/ (coef_mle_sips[3] - conf_val * std_error3) ,
    #              ## reverse 1/n
    #              coef_mle_sips[4] + conf_val * std_error4
    #              )
    #   lower <- c(coef_mle_sips[1] - conf_val * std_error1,
    #              coef_mle_sips[2] - conf_val * std_error2,
    #             exp(coef_mle_sips[2] - conf_val * std_error2),
    #              # exp(coef_mle[2]) - conf_val * std_error2_delta_method,
    #              coef_mle_sips[3] - conf_val * std_error2,
    #              # coef_mle_sips[3] - conf_val * std_error3,
    #              # 1/(coef_mle_sips[3] + conf_val * std_error3 ),
    #              ## reverse 1/n
    #              coef_mle_sips[4] - conf_val * std_error4
    #              )
    # 
    # #  table_result <- data.frame(Par = c("Qmax", "K", "n","sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
    #  #                            Upper=format(upper, scientific = T), Lower=format(lower, scientific = T))
    #   
    # 
    #   table_result <- data.frame(Par = c("Qmax", "Kprime", "K", "n","sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
    #                              Upper=format(upper, scientific = T), Lower=format(lower, scientific = T))
    #   
    #   # table_result <- data.frame(Par = c("Qmax", "Kprime", "K", "nprime","n","sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
    #   #                            Upper=format(upper, scientific = T), Lower=format(lower, scientific = T))
    #   # 
    # 
    # },
    # bordered=T)
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # output$hessian_num_original_inverse_sips <- renderTable({
    #   value <- inputVar()
    # 
    #   kth <- as.numeric(value[1])
    # 
    #   data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
    #   data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
    # 
    #   n_data_kth <- dim(data_kth)[1]
    #   x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
    #   initial_values <- data_meta[ kth, c(9,10)]
    #   location_error_vec <- NA
    #   colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
    # 
    # 
    #   ### 0. Commerical product
    # 
    # 
    #   ### 2. Freundlich model
    #   ## 2-1 initial value
    #   Y_trans <-  (1 / data_kth[-1,2] )
    #   X_trans <- (1 / data_kth[-1,1] )
    # 
    #   data_trans <- data.frame(q.mol.kg = Y_trans,
    #                            Cw.mol.l = X_trans)
    # 
    #   fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
    #   coef_lm <- fit_lm$coefficients
    #   Q_est_lm <- 1 / coef_lm[1]
    #   K_est_lm <- ( coef_lm[1] / coef_lm[2] )
    # 
    #   # Compute mse
    #   fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
    # 
    #   res_lm <- data_kth[-1,2] - 1/fitted_values_lm
    #   mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
    # 
    # 
    #   ### 2. MLE
    #   # 2-1 initial values from the previous linear model
    #   sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
    #   initial_value_mle <- c( 1/coef_lm[1], 
    #                           log( coef_lm[1] / coef_lm[2],),
    #                           1,
    #                           sqrt(sigma2_initial))
    # 
    #   # 3-2 fit a new model of reparameterization for Kd
    #   if(value[4] == "optim"){
    #     fit_mle_sips <- optim(fn = negLogLikelihood_sips,
    #                           par = initial_value_mle,
    #                           hessian = T,
    #                           method = value[5],
    #                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    #     coef_mle_sips <- fit_mle_sips$par
    # 
    #   } else if(value[4] == "nlm"){
    #     fit_mle_sips <- nlm(f = negLogLikelihood_sips,
    #                         p = initial_value_mle,
    #                         hessian = T,
    #                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
    #     )
    # 
    #     coef_mle_sips <- fit_mle_sips$estimate
    # 
    #   }
    # 
    # 
    # 
    #   solve(fit_mle_sips$hessian)
    # 
    # },
    # colnames = F,
    # digits=-2,
    # bordered = T
    # )



    
    output$table <- renderDataTable({
        value <- inputVar()

        kth <- as.numeric(value[1])

        data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
        data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]


        datatable(data_kth)
    })
} # end of server

# Run the application 
shinyApp(ui = ui, server = server)
