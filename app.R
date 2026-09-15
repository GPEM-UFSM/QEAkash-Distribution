# Aplicativo Shiny — Quantile Exponentiated Akash
# Dependências: install.packages(c("shiny","ggplot2","purrr","viridis","bslib","dplyr","plotly"))
# Execute no RStudio com Run App, ou: shiny::runApp("app_EAKASH.R")
# Todas as funções do app de referência foram adaptadas; tau é sempre explícito.

library(shiny)
library(ggplot2)
library(purrr)
library(viridis)
library(bslib)
library(dplyr)
library(plotly)
# Quantile Exponentiated Akash: mu = Q(tau), sigma > 0, tau fixado.
# F(y) = G(y)^alpha, alpha = log(tau)/log(G(mu)).
# G e g são a CDF e a densidade Akash do arquivo GAMLSS_EAKASH.R.
# Arquivo autossuficiente: não executa os exemplos ou a família GAMLSS original.
check_eakash <- function(mu, sigma, tau) {
  if (any(!is.finite(c(mu, sigma, tau))) || any(mu <= 0) ||
      any(sigma <= 0) || any(tau <= 0 | tau >= 1))
    stop("Use mu > 0, sigma > 0 e 0 < tau < 1, todos finitos.")
}
log1mexp <- function(a) {
  out <- a
  k <- !is.na(a) & a < -log(2)
  out[k] <- log1p(-exp(a[k]))
  out[!k] <- log(-expm1(a[!k]))
  out
}
logadd <- function(a, b) {
  m <- pmax(a,b)
  z <- m + log(exp(a-m) + exp(b-m))
  z[is.infinite(m)] <- m[is.infinite(m)]
  z
}
# A Akash é uma mistura de Gamma(1, sigma) e Gamma(3, sigma).
# Usar ambas as caudas evita a subtração instável 1 - G(y).
logG_akash <- function(y, sigma) {
  den <- logadd(2*log(sigma), log(2))
  a <- 2*log(sigma)-den
  b <- log(2)-den
  lo <- logadd(a+pgamma(y,1,rate=sigma,log.p=TRUE),
               b+pgamma(y,3,rate=sigma,log.p=TRUE))
  hi <- logadd(a+pgamma(y,1,rate=sigma,lower.tail=FALSE,log.p=TRUE),
               b+pgamma(y,3,rate=sigma,lower.tail=FALSE,log.p=TRUE))
  k <- !is.na(lo) & lo > -log(2)
  lo[k] <- log1mexp(hi[k])
  lo
}
alpha_eakash <- function(mu, sigma, tau) {
  check_eakash(mu,sigma,tau)
  a <- log(tau)/logG_akash(mu,sigma)
  if (any(!is.finite(a) | a <= 0))
    stop("Parâmetros fora da precisão numérica; reduza mu ou sigma.")
  a
}
dEAKASH <- function(y, mu, sigma, tau=0.5, log=FALSE) {
  a <- alpha_eakash(mu,sigma,tau)
  n <- max(length(y),length(a),length(sigma))
  y <- rep_len(y,n); a <- rep_len(a,n); sigma <- rep_len(sigma,n)
  ans <- rep(-Inf,n); ans[is.na(y)] <- NA_real_
  k <- which(is.finite(y) & y > 0)
  z <- y[k]; s <- sigma[k]
  lg <- 3*log(s)-logadd(2*log(s),log(2)) +
    logadd(0,2*log(z))-s*z
  ans[k] <- log(a[k])+lg+(a[k]-1)*logG_akash(z,s)
  k <- which(!is.na(y) & y == 0)
  ans[k] <- ifelse(a[k]<1,Inf,ifelse(a[k]>1,-Inf,
                   3*log(sigma[k])-logadd(2*log(sigma[k]),log(2))))
  if(log) ans else exp(ans)
}
pEAKASH <- function(q, mu, sigma, tau=0.5, lower.tail=TRUE, log.p=FALSE) {
  a <- alpha_eakash(mu,sigma,tau)
  lp <- a*logG_akash(q,sigma)
  if(!lower.tail) lp <- log1mexp(lp)
  if(log.p) lp else exp(lp)
}
sEAKASH <- function(q, mu, sigma, tau=0.5, log.p=FALSE)
  pEAKASH(q,mu,sigma,tau,lower.tail=FALSE,log.p=log.p)
hEAKASH <- function(y, mu, sigma, tau=0.5)
  exp(dEAKASH(y,mu,sigma,tau,log=TRUE)-sEAKASH(y,mu,sigma,tau,log.p=TRUE))
qEAKASH <- function(p, mu, sigma, tau=0.5, lower.tail=TRUE, log.p=FALSE) {
  a <- alpha_eakash(mu,sigma,tau)
  lp <- if(log.p) p else log(p)
  if(anyNA(lp) || any(lp>0)) stop("Probabilidades inválidas.")
  if(!lower.tail) lp <- log1mexp(lp)
  n <- max(length(lp),length(a),length(mu),length(sigma))
  lp <- rep_len(lp,n); a <- rep_len(a,n)
  mu <- rep_len(mu,n); sigma <- rep_len(sigma,n)
  vapply(seq_len(n),function(i) {
    if(lp[i]==-Inf) return(0)
    if(lp[i]==0) return(Inf)
    target <- lp[i]/a[i]
    # Raiz na escala log(y), sem intervalo superior fixo para y.
    f <- function(t) logG_akash(exp(t),sigma[i])-target
    center <- log(mu[i]); lo <- center; hi <- center
    while(f(lo)>0 && lo > -700) lo <- lo-4
    while(f(hi)<0 && hi < 700) hi <- hi+4
    exp(uniroot(f,c(lo-1,hi+1),tol=1e-10)$root)
  },numeric(1))
}
qEAKASH_numeric <- qEAKASH
rEAKASH <- function(n, mu, sigma, tau=0.5) {
  if(length(n)!=1 || !is.finite(n) || n<0 || n!=floor(n))
    stop("n deve ser um inteiro não negativo.")
  if(n==0) return(numeric())
  qEAKASH(runif(n),mu,sigma,tau)
}
reakash <- rEAKASH
bowley <- function(mu,sigma,tau=0.5) {
  q <- qEAKASH(c(.25,.5,.75),mu,sigma,tau)
  (q[3]-2*q[2]+q[1])/(q[3]-q[1])
}
moors <- function(mu,sigma,tau=0.5) {
  q <- qEAKASH(c(1,2,3,5,6,7)/8,mu,sigma,tau)
  (q[6]-q[4]+q[3]-q[1])/(q[5]-q[2])
}
loglik_eakash <- function(y,mu,sigma,tau=0.5) {
  if(any(!is.finite(y) | y<=0)) stop("A amostra deve ser positiva e finita.")
  sum(dEAKASH(y,mu,sigma,tau,log=TRUE))
}
# Busca multistart em log(mu), log(sigma); máximo numérico, não prova global.
fit_eakash <- function(y,tau,mu0,sigma0) {
  objective <- function(par) {
    val <- tryCatch(-loglik_eakash(y,exp(par[1]),exp(par[2]),tau),
                    error=function(e) Inf)
    if(is.finite(val)) val else 1e100
  }
  start <- expand.grid(mu=c(mu0,as.numeric(quantile(y,tau))),
                       sigma=sigma0*c(.3,1,3))
  fits <- lapply(seq_len(nrow(start)),function(i)
    tryCatch(optim(log(as.numeric(start[i,])),objective,method="Nelder-Mead",
                   control=list(maxit=2000,reltol=1e-10)),error=function(e) NULL))
  fits <- Filter(function(z) !is.null(z) && is.finite(z$value) && z$value<1e99, fits)
  if(!length(fits)) stop("A otimização falhou. Tente outra amostra ou parâmetros.")
  converged <- Filter(function(z) z$convergence==0,fits)
  if(length(converged)) fits <- converged
  best <- fits[[which.min(vapply(fits,function(z) z$value,numeric(1)))]]
  list(mu=exp(best$par[1]),sigma=exp(best$par[2]),loglik=-best$value,
       convergence=best$convergence)
}

# ============================================================================
# INTERFACE SHINY
# ============================================================================

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  withMathJax(),
  titlePanel("Quantile Exponentiated Akash Distribution — EAKASH"),
  h5("Assimetria, curtose, distribuição e log-verossimilhança"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tab_ativa == 'Log-verossimilhança'",
        h4("Simulação e estimação"),
        sliderInput("tau_l", "τ (fixo)", .01, .99, .5, step=.01),
        numericInput("mu_l", "μ verdadeiro = Q(τ)", 1, min=.01, max=20, step=.1),
        numericInput("sigma_l", "σ verdadeiro", 1, min=.01, max=10, step=.1),
        numericInput("n_l", "Tamanho da amostra", 200, min=10, max=5000, step=10),
        numericInput("seed_l", "Semente", 123, min=0, max=2147483646, step=1),
        actionButton("simulate", "Gerar amostra e estimar", class="btn-primary"),
        hr(),
        sliderInput("width_l", "Amplitude ao redor dos parâmetros (fator)",
                    1.2, 4, 2, step=.1),
        sliderInput("grid_l", "Resolução da superfície", 25, 80, 45, step=5),
        helpText("A amostra só muda ao clicar no botão. Alterar a amplitude ou a resolução mantém a mesma amostra."),
        downloadButton("sample_csv", "Baixar amostra")
      ),
      conditionalPanel(
        condition = "input.tab_ativa == 'Assimetria e curtose'",
        h4("Assimetria e curtose"),
        helpText("$$\\tau$$"),
        sliderInput("tau_ac", NULL, min = 0.01, max = 0.99, value = 0.5, step = 0.01),
        helpText("$$\\mu_{\\tau}$$"),
        sliderInput("mu_max", NULL, min = 1, max = 10, value = 5, step = 0.1),
        helpText("$$\\sigma_{max}$$"),
        sliderInput("sigma_max", NULL, min = 1, max = 8, value = 5, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.tab_ativa == 'Distribuição'",
        h4("Distribuição"),
        helpText("$$\\tau$$"),
        sliderInput("tau_d", NULL, min = 0.01, max = 0.99, value = 0.5, step = 0.01),
        helpText("$$\\mu_{\\tau}$$"),
        sliderInput("mu_d", NULL, min = 0.5, max = 10, value = 1, step = 0.1),
        helpText("$$\\sigma$$"),
        sliderInput("sigma_d", NULL, min = 0.1, max = 5, value = 1, step = 0.1),
        helpText("$$y_{max}$$"),
        sliderInput("y_max", NULL, min = 1, max = 10, value = 5, step = 0.1)
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "tab_ativa", selected = "Distribuição",
        tabPanel(
          "Distribuição",
          fluidRow(
            column(6, plotOutput("dist_dens")),
            column(6, plotOutput("dist_cdf")),
            column(6, plotOutput("dist_surv")),
            column(6, plotOutput("dist_haz"))
          )
        ),
        tabPanel(
          "Assimetria e curtose",
          plotOutput("bowley"),
          br(),
          plotOutput("moors")
        ),
        tabPanel("Log-verossimilhança",
          helpText("τ é fixado. μ é o quantil de ordem τ, não necessariamente a mediana. O ponto vermelho indica o máximo numérico encontrado; o losango azul indica os parâmetros geradores."),
          verbatimTextOutput("ll_status"),
          tableOutput("ll_estimates"),
          plotlyOutput("ll_surface", height="530px"),
          plotOutput("ll_contour", height="400px"),
          sliderInput("hist_bins", "Número de classes do histograma", 10, 80, 30, step=1),
          checkboxInput("hist_expand", "Expandir a escala vertical para incluir todo o pico avaliado", FALSE),
          plotOutput("ll_sample", height="380px"),
          helpText("A superfície representa ℓ(μ, σ) = Σ log f(yᵢ; μ, σ, τ). A busca usa múltiplos valores iniciais; convergência não garante máximo global.")
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output) {
  

  experiment <- eventReactive(input$simulate, {
    validate(need(is.finite(input$mu_l) && input$mu_l>0 && input$mu_l<=20,
                  "Use 0 < μ ≤ 20."),
             need(is.finite(input$sigma_l) && input$sigma_l>0 && input$sigma_l<=10,
                  "Use 0 < σ ≤ 10."),
             need(is.finite(input$n_l) && input$n_l>=10 && input$n_l<=5000 &&
                  input$n_l==floor(input$n_l), "Use n inteiro entre 10 e 5000."),
             need(is.finite(input$seed_l) && input$seed_l>=0 &&
                  input$seed_l<=2147483646 && input$seed_l==floor(input$seed_l),
                  "Use uma semente inteira válida."))
    withProgress(message="Simulando e maximizando a log-verossimilhança", value=.1, {
      tryCatch({
        set.seed(input$seed_l)
        y <- rEAKASH(input$n_l,input$mu_l,input$sigma_l,input$tau_l)
        incProgress(.3)
        fit <- fit_eakash(y,input$tau_l,input$mu_l,input$sigma_l)
        list(y=y,fit=fit,mu=input$mu_l,sigma=input$sigma_l,tau=input$tau_l,
             seed=input$seed_l)
      }, error=function(e) { validate(need(FALSE,conditionMessage(e))) })
    })
  }, ignoreNULL=FALSE)
  llgrid <- reactive({
    e <- experiment(); w <- input$width_l
    mus <- exp(seq(log(min(e$mu,e$fit$mu)/w),log(max(e$mu,e$fit$mu)*w),
                   length.out=input$grid_l))
    sigmas <- exp(seq(log(min(e$sigma,e$fit$sigma)/w),log(max(e$sigma,e$fit$sigma)*w),
                      length.out=input$grid_l))
    # plotly espera linhas em sigma e colunas em mu.
    z <- outer(sigmas,mus,Vectorize(function(s,m)
      tryCatch(loglik_eakash(e$y,m,s,e$tau),error=function(err) NA_real_)))
    z[!is.finite(z)] <- NA_real_
    list(mu=mus,sigma=sigmas,z=z)
  })
  output$ll_status <- renderText({
    e <- experiment(); g <- llgrid()
    text <- sprintf("Amostra: n = %d | semente = %d | τ = %.2f\nConvergência: %s | ℓ máximo encontrado = %.6f",
      length(e$y),e$seed,e$tau,if(e$fit$convergence==0) "sim" else "não",e$fit$loglik)
    if(any(g$z > e$fit$loglik+1e-5,na.rm=TRUE))
      text <- paste(text,"Atenção: a grade contém valor superior ao ajuste; o máximo não foi confirmado.",sep="\n")
    if(anyNA(g$z)) text <- paste(text,"Alguns pontos da superfície estão fora da precisão numérica e foram omitidos.",sep="\n")
    text
  })
  output$ll_estimates <- renderTable({
    e <- experiment()
    data.frame(Parâmetro=c("μ","σ"),Verdadeiro=c(e$mu,e$sigma),
               Estimado=c(e$fit$mu,e$fit$sigma),check.names=FALSE)
  },digits=6)
  output$ll_surface <- renderPlotly({
    e <- experiment(); g <- llgrid()
    true_ll <- loglik_eakash(e$y,e$mu,e$sigma,e$tau)
    fig <- plot_ly(x=g$mu,y=g$sigma,z=g$z,type="surface",colors=viridis::viridis(256),
                   name="Log-verossimilhança",showlegend=FALSE,
                   colorbar=list(title="ℓ",x=1.02,y=.55,len=.8,thickness=16))
    fig <- add_trace(fig,x=e$fit$mu,y=e$fit$sigma,z=e$fit$loglik,
      type="scatter3d",mode="markers",inherit=FALSE,name="Máximo numérico",showlegend=TRUE,
      marker=list(color="red",size=6),
      hovertemplate="μ = %{x:.6f}<br>σ = %{y:.6f}<br>ℓ = %{z:.6f}<extra>Máximo numérico</extra>")
    fig <- add_trace(fig,x=e$mu,y=e$sigma,z=true_ll,type="scatter3d",
      mode="markers",inherit=FALSE,name="Parâmetros verdadeiros",showlegend=TRUE,
      marker=list(color="dodgerblue",size=5,symbol="diamond"))
    layout(fig,
      scene=list(xaxis=list(title="μ"),yaxis=list(title="σ"),
                 zaxis=list(title="ℓ(μ, σ)")),
      showlegend=TRUE,
      legend=list(orientation="h",x=.5,xanchor="center",y=-.12,yanchor="top",
                  bgcolor="rgba(255,255,255,0.9)",font=list(size=12)),
      margin=list(l=10,r=85,b=105,t=25))
  })
  output$ll_contour <- renderPlot({
    e <- experiment(); g <- llgrid()
    df <- expand.grid(sigma=g$sigma,mu=g$mu); df$ll <- as.vector(g$z)
    ggplot(df,aes(mu,sigma,z=ll)) + geom_contour_filled(bins=18,na.rm=TRUE) +
      annotate("point",x=e$fit$mu,y=e$fit$sigma,
               color="red",size=3) +
      annotate("point",x=e$mu,y=e$sigma,
               color="dodgerblue",shape=18,size=3) +
      labs(x=expression(mu),y=expression(sigma),fill="ℓ",title="Curvas de nível") +
      theme_minimal()
  })
  output$ll_sample <- renderPlot({
    e <- experiment()
    # Classes explícitas: cada altura é frequência / (n * largura da classe).
    # Toda a amostra participa do histograma; a área total é 1.
    xmax <- max(e$y)*1.05
    breaks <- seq(0,xmax,length.out=input$hist_bins+1L)
    hh <- hist(e$y,breaks=breaks,plot=FALSE,include.lowest=TRUE)
    bars <- data.frame(left=head(hh$breaks,-1),right=tail(hh$breaks,-1),
                       density=hh$density)
    # A densidade pode divergir em zero. Uma grade também logarítmica
    # resolve a curvatura perto da origem sem ligar só dois pontos distantes.
    xmin <- max(min(e$y)*.1,.Machine$double.xmin)
    x <- sort(unique(c(exp(seq(log(xmin),log(xmax),length.out=600)),
                       seq(xmin,xmax,length.out=800))))
    df <- data.frame(y=x,
                     ajustada=dEAKASH(x,e$fit$mu,e$fit$sigma,e$tau))
    df <- df[is.finite(df$ajustada),,drop=FALSE]
    alpha <- alpha_eakash(e$fit$mu,e$fit$sigma,e$tau)
    # Zoom por coordenadas: não descarta dados, não renormaliza a densidade.
    ymax <- 1.25*max(hh$density)
    clipped <- !isTRUE(input$hist_expand) && any(df$ajustada>ymax)
    note <- if(alpha<1) {
      sprintf("α = %.4f < 1: a densidade ajustada tende a infinito quando y se aproxima de zero.",alpha)
    } else {
      sprintf("α = %.4f. Histograma normalizado: área total igual a 1.",alpha)
    }
    if(clipped) note <- paste(note,
      "Escala ampliada para o histograma: a parte da curva acima do painel foi recortada visualmente.",sep="\n")
    if(isTRUE(input$hist_expand) && alpha<1) note <- paste(note,
      "O pico exibido depende do menor y positivo avaliado; o limite em zero é infinito.",sep="\n")
    fig <- ggplot() +
      geom_rect(data=bars,aes(xmin=left,xmax=right,ymin=0,ymax=density),
                fill="#B9CFDF",color="white",linewidth=.3) +
      geom_line(data=df,aes(x=y,y=ajustada,color="Densidade ajustada"),linewidth=1) +
      scale_color_manual(values=c("Densidade ajustada"="#D62728")) +
      labs(x="y",y="Densidade",color=NULL,
           title="Histograma da amostra simulada e densidade ajustada",
           subtitle=sprintf("Estimativas: μ = %.4f | σ = %.4f | τ fixo = %.2f",
                            e$fit$mu,e$fit$sigma,e$tau),caption=note) +
      theme_minimal(base_size=12) +
      theme(legend.position="top",plot.caption=element_text(hjust=0))
    if(isTRUE(input$hist_expand)) {
      fig + coord_cartesian(xlim=c(0,xmax),expand=FALSE)
    } else {
      fig + coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax),expand=FALSE)
    }
  })
  output$sample_csv <- downloadHandler(filename=function() "amostra_EAKASH.csv",
    content=function(file) {
      e <- experiment()
      write.csv(data.frame(y=e$y,mu=e$mu,sigma=e$sigma,tau=e$tau,seed=e$seed),
                file,row.names=FALSE)
    })

  # Dados para Bowley - varia mu
  db <- reactive({
    req(input$mu_max, input$tau_ac, input$sigma_max)
    mu_seq <- seq(0.1, input$mu_max, length.out = 10)
    sigma_seq <- seq(0.1, input$sigma_max, length.out = 50) # Reduzido levemente para ganho de performance
    
    map_df(mu_seq, function(mu) {
      bowley_vals <- sapply(sigma_seq, function(s) {
        tryCatch(bowley(mu = mu, sigma = s, tau = input$tau_ac),
                 error = function(e) NA_real_)
      })
      data.frame(
        sigma = sigma_seq,
        v = bowley_vals,
        mu_label = factor(sprintf("%.2f", mu), levels = sprintf("%.2f", mu_seq))
      )
    })
  })
  
  # Dados para Moors - varia mu
  dm <- reactive({
    req(input$mu_max, input$tau_ac, input$sigma_max)
    mu_seq <- seq(0.1, input$mu_max, length.out = 10)
    sigma_seq <- seq(0.1, input$sigma_max, length.out = 50)
    
    map_df(mu_seq, function(mu) {
      moors_vals <- sapply(sigma_seq, function(s) {
        tryCatch(moors(mu = mu, sigma = s, tau = input$tau_ac),
                 error = function(e) NA_real_)
      })
      data.frame(
        sigma = sigma_seq,
        v = moors_vals,
        mu_label = factor(sprintf("%.2f", mu), levels = sprintf("%.2f", mu_seq))
      )
    })
  })
  
  # Dados para Distribuição
  dd <- reactive({
    req(input$mu_d, input$sigma_d, input$tau_d, input$y_max)
    y_seq <- seq(0.001, input$y_max, length.out = 100)
    
    df <- data.frame(
      y = y_seq,
      dens = dEAKASH(y_seq, input$mu_d, input$sigma_d, input$tau_d),
      cdf = pEAKASH(y_seq, input$mu_d, input$sigma_d, input$tau_d),
      surv = sEAKASH(y_seq, input$mu_d, input$sigma_d, input$tau_d),
      haz = hEAKASH(y_seq, input$mu_d, input$sigma_d, input$tau_d)
    )
    df[is.finite(df$dens) & is.finite(df$cdf) & is.finite(df$surv) & is.finite(df$haz), , drop = FALSE]
  })
  
  # Função auxiliar para plotar Bowley/Moors
  plot_gg <- function(df, ylab) {
    req(df, nrow(df) > 0)
    fin <- df[is.finite(df$v), ]
    req(nrow(fin) > 0)
    
    n_colors <- length(unique(fin$mu_label))
    
    ggplot(fin, aes(x = sigma, y = v, colour = mu_label)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(
        values = viridis::mako(n_colors, alpha = 0.8, begin = 0.1, end = 0.9),
        name = expression(mu[tau])
      ) +
      labs(y = ylab, x = expression(sigma), title = "") +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box = "vertical",
        legend.position = "right"
      ) +
      guides(colour = guide_legend(
        direction = "vertical", title.position = "top",
        override.aes = list(linewidth = 1)
      )) +
      coord_cartesian(ylim = c(min(fin$v, na.rm = TRUE), max(fin$v, na.rm = TRUE)), xlim = c(0, input$sigma_max))
  }
  
  output$bowley <- renderPlot(plot_gg(db(), "Bowley Skewness"))
  output$moors <- renderPlot(plot_gg(dm(), "Moors Kurtosis"))
  
  # Gráficos de Distribuição
  output$dist_dens <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = dens)) +
      geom_line(color = "#440154", linewidth = 1) +
      labs(title = "Probability Density Function", x = "y", y = "f(y)") +
      theme_minimal()
  })
  
  output$dist_cdf <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = cdf)) +
      geom_line(color = "#3B528B", linewidth = 1) +
      labs(title = "Cumulative Distribution Function", x = "y", y = "F(y)") +
      theme_minimal()
  })
  
  output$dist_surv <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = surv)) +
      geom_line(color = "#21918C", linewidth = 1) +
      labs(title = "Survival Function", x = "y", y = "S(y)") +
      theme_minimal()
  })
  
  output$dist_haz <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = haz)) +
      geom_line(color = "brown", linewidth = 1) +
      labs(title = "Hazard Function", x = "y", y = "h(y)") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)

