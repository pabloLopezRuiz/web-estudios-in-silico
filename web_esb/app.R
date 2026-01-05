library(shiny)
library(bslib)
# library(DT)
# library(reactable)
library(ggplot2)
library(plotly)
# library(heatmaply)
# library(limma)
# library(visNetwork)
# library(DESeq2)
# library(mixOmics)
# library(STRINGdb)

 
load("GSE124799/gse.RData")
load("data/coldata.RData")
load("data/dmtx.RData")

# dea results
load("GSE124799/DE_Sed_vs_Ex.RData")
load("GSE124799/DE_Sed_vs_ExSed.RData")
load("GSE124799/DE_Ex_vs_ExSed.RData")
load("data/corr_matrix.RData")

# PCA files
load("data/pca_res.RData")
load("data/mixomics_pca.RData")

# func prof
load("data/sve_fun.RData")
load("data/eves_fun.RData")
load("data/sves_fun.RData")

# UI
ui <- page_navbar(
  title = "Análisis de RNA-Seq: GSE124799",
  theme = bs_theme(
    bootswatch = "minty",
    base_font = font_google("Inter"),
    "card-cap-bg" = "#f8f9fa"
  ),

  # CSS Personalizado para márgenes globales
  header = tags$head(
    tags$style(HTML("
      body {
        background-color: #f4f7f6 !important;
      }
      .container-fluid {
        padding: 2rem 5% !important;
      }
      .card {
        box-shadow: 0 10px 25px rgba(0,0,0,0.05);
        border: none !important;
        margin-bottom: 24px;
        transition: transform 0.2s ease-in-out;
      }
      .card:hover {
        transform: translateY(-5px);
        box-shadow: 0 15px 30px rgba(0,0,0,0.1);
      }
      .card-header {
        background-color: #ffffff !important;
        border-bottom: 1px solid #f0f0f0 !important;
        font-weight: 600;
        letter-spacing: 0.5px;
      }
      .navbar {
        box-shadow: 0 2px 10px rgba(0,0,0,0.1);
      }
    "))
  ),

  # Tab 1: Overview
  nav_panel(
    title = "Resumen",
    div(
      h4("Estudio Elegido"),
      p(
        "Análisis de RNA-seq en ratones C57BL/6 de 7-8 semanas con acceso a ruedas de ejercicio",
        " durante 6 semanas (Ex), individuos sedentarios (Sed) y un grupo con ejercicio",
        " seguido de 3 semanas de inactividad (ExSed). Perfiles LSK por Illumina NextSeq."
      ),
      p(strong("Modelo experimental: ratones C57BL/6")),
      p(strong("Muestras: 9")),
      p(strong("Genes:  34752")),
      hr(),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Sample distribution: Sex"),
          plotlyOutput("hist_sex")
        ),
        card(
          card_header("Sample distribution: Experimental group"),
          plotlyOutput("hist_group")
        ),
        card(
          card_header("Density distribution"),
          plotOutput("plot_density")
        ),
        card(
          card_header("Boxplot"),
          input_switch(id = "norm", label = "Datos normalizados", value = F),
          plotlyOutput("box_prenorm")
        ),
      ),
      card(
        card_header("Correlation"),
        plotlyOutput("corr_plot")
      ),
      accordion(
        open = FALSE,
        class = "mb-4",
        accordion_panel(
          title = "Observando los datos...",
          icon = icon("book-open"),
          p(
            "En un primer vistazo a los datos, observamos que la distribución",
            " no es equitativa en cuanto a sexo, ya que el estudio cuenta",
            " únicamente con ratones macho. No obstante, el número de individuos",
            " de cada grupo experimental sí que es equitativo (Sed, Ex y ExSed)."
          ),
          p(
            "Respecto a la densidad de conteos y la distribución, los datos presentan el",
            " perfil característico de los experimentos de RNA-seq: muchos genes",
            " con baja expresión o conteos nulos, y un grupo reducido de genes altamente",
            " expresados que sesgan la distribución inicial. La aplicación de la",
            " la transformación VST en los boxplots permite mitigar el efecto",
            " de los 'outliers'"
          ),
          p(
            "Finalmente, el análisis de correlación se observa una homogeneidad global",
            " muy elevada (coeficientes entre 0.85 y 1), lo que confirma que la identidad",
            " transcriptómica de las células se mantiene estable pese a estudiar",
            " condiciones distintas. Se observa un cluster entre las muestras GSM3554211,",
            " GSM3554210 (ExSed) y GSM3554204 (Sed). Este agrupamiento sugiere que",
            " el perfil inducido por el ejercicio parece desvanecerse tras el",
            " periodo de inactividad. Vemos también que la GSM3554209 (ExSed)",
            " se encuentra relacionada con las anteriores además de con la",
            " GSM3554203 y la GSM3554206, lo cual parece dar mensajes",
            " contradictorios porque esta última pertenece al grupo Ex. Puede",
            " que estas diferencias observadas provengan de ruido técnico que",
            " no se ha podido filtar."
          )
        )
      )
    ),
  ),


  # --- TAB: DIM RED ---

  nav_panel(
    title = "PCA",
    div(
      card(
        card_header("PCA plot"),
        plotOutput("pca_plot")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Varianza explicada"),
          plotOutput("pca_var")
        ),
        card(
          card_header("Weights plot"),
          layout_columns(
            card(
              card_header("Número de genes mostrado"),
              sliderInput("ngenes",
                "",
                min = 5,
                max = 1000,
                value = 50,
                tick = F,
              ),
            ),
            card(
              card_header("Componente principal"),
              sliderInput("comp",
                "",
                min = 1,
                max = 5,
                value = 1,
                tick = F,
              ),
            ),
          ),
          card(plotOutput("pca_wei"))
        ),
      ),
      accordion(
        open = FALSE,
        class = "mb-4",
        accordion_panel(
          title = "Observando la PCA...",
          icon = icon("book-open"),
          p(
            "En la proyección de las componentes 1 y 2, se observa una",
            " agrupación vaga entre Sed y ExSed, diferenciados principalmente",
            " por PC1. El grupo Ex presenta mayor dispersión, sin formar un",
            " cluster definido. Es inesperado que el grupo ExSed no se sitúe",
            " en un punto medio entre Sed y Ex"
          ),
          p(
            "La componente PC1 captura casi el 40% de la varianza,",
            " predominando sobre el resto. Esta concentración de la",
            " varianza sugiere que la variabilidad depende de pocos factores.",
            " Dado que el PC1 discrimina entre ExSed y Sed de Ex, es probable",
            " que represente una señal biológica vinculada al ejercicio."
          ),
          p(
            "En los pesos del PC1, destaca un conjunto de 10 a 20",
            " genes con alta contribución, indicando que la fuente",
            " principal está determinada por pocos transcritos. En la",
            " componente 2, hay una polarización global, todos",
            " los genes tienen pesos negativos, excepto",
            " ENSMUSG00000066632 y ENSMUSG00000095079. Esto indica",
            " que estos dos genes actúan de forma antagónica a la tendencia",
            " general."
          )
        )
      )
    ),
  ),

  # --- TAB: DEA ---
  nav_panel(
    title = "DEA",
    div(
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Código usado"),
          p(
            "Para el análisis de expresión diferencial hemos usado el paquete",
            " DESeq2 con el siguiente código"
          ),
          markdown(
            "```R
            # Cargo las librerías
            library(DESeq2)

            # Cargo el dataset (en mi github)
            load(\"GSE124799/gse.RData\")

            # Elimino los genes con pocos conteos
            keep <- rowSums(assay(g) >= 2) >= 3
            g_sel <- g[keep, ]

            # Establezco el nivel sedentario como control
            g_sel$group <- relevel(g_sel$group, ref = \"Sed\")

            # Transformo mis datos a DESeq, los normalizo y hago el modelo
            ddsSE <- DESeqDataSet(g_sel, design = ~group)
            dsNorm <- estimateSizeFactors(ddsSE, type=\"ratio\")
            dds <- DESeq(dsNorm)

            # Guardo los distintos contrastes
            SvE <- results(dds, contrast = c(\"group\", \"Sed\", \"Ex\"))
            SvEs <- results(dds, contrast = c(\"group\", \"Sed\", \"ExSed\"))
            EvEs <- results(dds, contrast = c(\"group\", \"Ex\", \"ExSed\"))
            ```"
          ),
        ),
        card(
          style = "border: none;",
          radioButtons(
            "ctr",
            "Contraste",
            choices = c(
              "Sed vs Ex"    = "sve",
              "Ex vs ExSed"  = "eves",
              "Sed vs ExSed" = "sves"
            ),
            selected = "sve"
          ),
          card(
            card_header("Matriz de diseño experimental"),
            reactable::reactableOutput("table")
          ),
          uiOutput("form")
        ),
      ),
      card(
        card_header("Top genes"),
        layout_sidebar(
          sidebar = sidebar(
            sliderInput(
              "ptopgenes",
              "P valor mínimo",
              min = 0.001,
              max = 0.5,
              value = 0.05,
              ticks = F
            ),
            sliderInput(
              "ntopgenes",
              "Número de genes mostrado",
              min = 5,
              max = 222,
              value = 20,
              ticks = F
            ),
            selectInput(
              "otopgenes",
              "Selecciona el tipo de ordenación:",
              list(
                "Ascendente" = "asc",
                "Descentente" = "des",
                "valor absoluto" = "abs"
              )
            ),
          ),
          plotlyOutput("topgenes")
        ),
      ),
      card(
        card_header("Volcano Plot"),
        layout_sidebar(
          sidebar = sidebar(
            sliderInput(
              "pvlc",
              "P valor mínimo",
              min = 0.001,
              max = 0.5,
              value = 0.05,
              ticks = F
            ),
            sliderInput(
              "lvlc",
              "LogFoldChange mínimo",
              min = 0,
              max = 25,
              value = 2,
              ticks = F
            ),
          ),
          plotlyOutput("vlc_plot")
        ),
      ),
      accordion(
        open = FALSE,
        class = "mb-4",
        accordion_panel(
          title = "Observando el DEA...",
          icon = icon("book-open"),
          p(
            "En línea con los análisis previos, los resultados del DEA",
            " confirman que las mayores divergencias transcriptómicas se",
            " producen al comparar el ejercicio activo (Ex) frente a los",
            " estados de inactividad (Sed o ExSed)."
          ),
          p(
            "En el primer contraste (Ex ~ ExSed), bajo un umbral de",
            " significancia de 0.05, identificamos 34 genes diferencialmente",
            " expresados con logFC de hasta 20. Este resultado evidencia",
            " que el perfil transcripcional de los individuos activos es",
            " significativamente distinto al de aquellos que han dejado de",
            " ejercitarse, lo que implica un cambio del transcriptoma tras",
            " el cese del mismo"
          ),
          p(
            "Por el contrario, la comparación Sed ~ ExSed revela solo 8",
            " genes diferencialmente expresados. De este grupo reducido,",
            " únicamente ENSMUSG00000094103 presenta un logFC superior a 20.",
            " Es interesante apuntar que este gen se encuentra entre los 30",
            " principales contribuyentes de la PC1, aunque no forma parte",
            " del grupo reducido de genes que contribuyen mayoritariamente",
            " a la variabilidad de esta componente."
          ),
          p(
            "Como conclusión, estos contrastes ratifican que el ejercicio",
            " físico provoca alteraciones en el transcriptoma de las",
            " células LSK de ratón. Sin embargo, la mínima diferencia",
            " observada entre los grupos Sed y ExSed sugiere que estos",
            " efectos son desaparecen en gran medida tras un periodo de",
            " inactividad, devolviendo al sistema a un estado molecular",
            " muy cercano al sedentario."
          )
        )
      )
    ),
  ),

  # --- TAB: CARACT FUNC ---
  nav_panel(
    title = "Caracterización funcional",
    div(
      layout_columns(
        col_widths = c(2, 10),
        card(
          card_header("Contraste:"),
          radioButtons(
            "ctr_fun",
            "Contraste",
            choices = c(
              "Sed vs Ex"    = "sve_fun",
              "Ex vs ExSed"  = "eves_fun",
              "Sed vs ExSed" = "sves_fun"
            ),
            selected = "sve_fun"
          ),
        ),
        card(
          card_header("String"),
          visNetwork::visNetworkOutput("string")
        ),
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Barplot Top genes"),
          plotOutput("ora_bar_top")
        ),
        card(
          card_header("Barplot Bottom genes"),
          plotOutput("ora_bar_bot")
        ),
        card(
          card_header("GSEA GO"),
          plotOutput("gsea_go")
        ),
        card(
          card_header("GSEA KEGG"),
          plotOutput("gsea_kegg")
        ),
      ),
      accordion(
        open = FALSE,
        class = "mb-4",
        accordion_panel(
          title = "Observando el análisis funcional...",
          icon = icon("book-open"),
          p(
            "El análisis de enriquecimiento funcional revela que los",
            " genes sobreexpresados en respuesta al ejercicio (Sed ~ Ex)",
            " están vinculados a rutas de señalización y supervivencia,",
            " destacando la regulación de la apoptosis y la señalización",
            " por factores de crecimiento como EGF y ERBB."
          ),
          p(
            "En el contraste Ex ~ ExSed, los genes sobreexpresados en el",
            " grupo de retirada muestran una asociación con la respuesta",
            " inmune (linfocitos B, granulocitos y neutrófilos). Esto",
            " sugiere un efecto rebote donde el cese del ejercicio activa",
            " programas proinflamatorios que el deporte mantenía",
            " reprimidos, o bien que el sedentarismo induce estos efectos",
            " de forma acelerada."
          ),
          p(
            "En la comparación Sed ~ ExSed, se observa de nuevo la",
            " activación de rutas inmunitarias como el estallido",
            " respiratorio. Este hallazgo apoya la hipótesis de que es la",
            " transición a la inactividad la que provoca esta respuesta",
            " inmune, más que el estado sedentario en sí mismo."
          ),
          p(
            "El análisis GSEA refuerza estas tendencias. En Sed ~ Ex,",
            " destaca la activación de la vía 'Hippo', conocida por",
            " inhibir la proliferación celular. En Ex ~ ExSed, vemos un",
            " aumento en apoptosis, fosforilación oxidativa y actividad",
            " ribosomal. El cese del ejercicio eleva la señalización por",
            " TNF y el metabolismo general, lo que podría vincularse a la",
            " inactivación del freno que suponía la vía Hippo. Para el",
            " último contraste vemos que tanto en GO como en KEGG destacan",
            " las rutas relacionadas con el aumento del metabolismo como la",
            " síntesis de ATP o el aumento de la fosforilación oxidativa, de",
            " nuevo apoyando la hipótesis presentada anteriormente"
          ),
          p(
            "Finalmente, la red de interacción de STRING se presenta en",
            " los 3 contrastes como una pequeña red bastante",
            " interconectada con 1 o 2 cluster de pocos genes, apoyando",
            " la idea inicial de que la mayoría de la variación proviene",
            " de unos pocos genes, que ahora podemos ver que se",
            " encuentran bastante interconectados."
          )
        )
      )
    ),
  ),
)


# Server
server <- function(input, output, session) {
  # Definicion del contraste elegido en cada momento
  df_contrast <- reactive({
    switch(input$ctr,
      sve  = SvE,
      eves = EvEs,
      sves = SvEs
    )
  })

  contrast_fun <- reactive({
    switch(input$ctr_fun,
      sve_fun  = sve_fun,
      eves_fun = eves_fun,
      sves_fun = sves_fun
    )
  })


  # --- Overview Tab ---
  output$hist_sex <- renderPlotly({
    p <- ggplot(coldata, aes(x = gender)) +
      scale_x_discrete(drop = FALSE) +
      geom_bar(width = 0.3, fill = "salmon") +
      xlab("Sexo") +
      ylab("Cantidad") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

    ggplotly(p)
  })

  output$hist_group <- renderPlotly({
    p <- ggplot(coldata, aes(x = group, fill = group)) +
      geom_bar(width = 0.3) +
      xlab("Grupo experimental") +
      ylab("Cantidad") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

    ggplotly(p)
  })

  output$plot_density <- renderPlot({
    limma::plotDensities(
      object = log2(SummarizedExperiment::assay(g) + 1),
      main = "Densidad de conteos (log2 + 1)",
      legend = F,
    )
    legend("topright", legend = colnames(g), col = 1:ncol(g), lty = 1)
  })

  output$box_prenorm <- renderPlotly({
    if (input$norm) {
      load("data/boxplot_prenorm.RData")
      ggplotly(bxpl)
    } else {
      load("data/boxplot_norm.RData")
      ggplotly(bxpl)
    }
  }) # -> OPTIMIZA

  output$corr_plot <- renderPlotly({
    # Correlation
    heatmaply::heatmaply(
      cor_matrix,
      ColSideColors = coldata$group,
      colors = colorRampPalette(c("navy", "white", "red"))(50),
      plot_method = "plotly"
    )
  })

  output$pca_plot <- renderPlot({
    mixOmics::plotIndiv(X,
      group = SummarizedExperiment::colData(g)$group,
      pch = "circle",
      legend = TRUE,
    )
  })

  output$pca_var <- renderPlot({
    PCAtools::screeplot(pca.res)
  })

  output$pca_wei <- renderPlot({
    mixOmics::plotLoadings(X, comp = input$comp, ndisplay = input$ngenes)
  })

  # --- DEA ---

  output$table <- reactable::renderReactable({
    con_code <- "
  function(rowInfo, state) {
    // Filas a colorear y sus colores
    const rowColors = {
      0: '#fae134',  // fila 1
      1: '#fae134',  // fila 2
      2: '#fae134',  // fila 3
      3: '#c2fa34',  // fila 4
      4: '#c2fa34',  // fila 5
      5: '#c2fa34',  // fila 6
      6: '#ff4d5b',  // fila 7
      7: '#ff4d5b',  // fila 8
      8: '#ff4d5b'   // fila 9
    };

    // Si el índice está definido en rowColors, retorna el color
    if (rowColors.hasOwnProperty(rowInfo.index)) {
      return { backgroundColor: rowColors[rowInfo.index] };
    } else {
      return { backgroundColor: '#ffffff' }; // resto de filas blanco
    }
  }
  "

    reactable::reactable(
      dmtx,
      rowStyle = reactable::JS(con_code)
    )
  })

  output$form <- renderUI({
    switch(input$ctr,
      sve = {
        top_color <- "#fae134"
        bottom_color <- "#c2fa34"
        top_name <- "Sed"
        bottom_name <- "Ex"
      },
      eves = {
        top_color <- "#c2fa34"
        bottom_color <- "#fa3444"
        top_name <- "Ex"
        bottom_name <- "ExSed"
      },
      sves = {
        top_color <- "#fae134"
        bottom_color <- "#fa3444"
        top_name <- "Sed"
        bottom_name <- "ExSed"
      }
    )
    div(
      style = "font-size: 24px; display: flex; align-items: center;",
      tags$span(
        style = paste0(
          "background-color:", top_color,
          "; color: black; padding:3px 8px; border-radius:4px; margin-right:5px;"
        ),
        paste0("\u2211 ", top_name)
      ),
      tags$span(
        style = "margin:0 20px;",
        "~"
      ),
      tags$span(
        style = paste0(
          "background-color:", bottom_color,
          "; color: black; padding:3px 8px; border-radius:4px; margin-right:5px;"
        ),
        paste0("\u2211 ", bottom_name)
      )
    )
  })

  output$topgenes <- renderPlotly({
    df <- df_contrast()
    df$genes <- rownames(df)

    # Seleccionamos por padj y preparamos el dataframe
    df <- as.data.frame(df[which(df$padj < input$ptopgenes), ])
    df <- df[, c("genes", "log2FoldChange", "padj")]

    # Arrange
    if (input$otopgenes == "asc") {
      df <- df[order(df$log2FoldChange), ]
    }
    if (input$otopgenes == "des") {
      df <- df[order(-df$log2FoldChange), ]
    }
    if (input$otopgenes == "abs") {
      df <- df[order(-abs(df$log2FoldChange)), ]
    }

    df$genes <- factor(df$genes, levels = df$genes)

    # Head
    df <- head(df, input$ntopgenes)

    # Plot
    fig <- ggplot(df, aes(x = genes, y = log2FoldChange)) +
      geom_col() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 7))

    ggplotly(fig)
  })

  output$vlc_plot <- renderPlotly({
    df <- df_contrast()
    sig_limit <- input$pvlc
    lfc_limit <- input$lvlc

    df[, "Significative"] <- rep("none", length(rownames(df)))

    df[which(df$padj < sig_limit &
      abs(df$log2FoldChange) < lfc_limit), "Significative"] <- "p.adjusted"
    df[which(df$padj < sig_limit &
      abs(df$log2FoldChange) > lfc_limit), "Significative"] <- "p.adjusted + logFC"

    ggplot(data = as.data.frame(df), aes(y = -log10(padj), x = log2FoldChange)) +
      geom_point(aes(color = Significative), na.rm = TRUE) +
      geom_vline(
        xintercept = c(-lfc_limit, lfc_limit),
        linetype = "dashed", color = "#E6007B"
      ) +
      geom_hline(
        yintercept = -log10(sig_limit),
        linetype = "dotted", color = "#734F77"
      ) +
      scale_color_manual(values = c("#A19D9F", "#E6007B", "#009E73")) +
      xlab("logFC") +
      ylab("-log10(p.adjusted)")
  })

  # --- GSA ---

  output$string <- visNetwork::renderVisNetwork({
    contrast_fun()$Network
  })

  output$ora_bar_top <- renderPlot({
    oratop <- contrast_fun()$ORA_Top
    barplot(oratop, showCategory = 5)
  })

  output$ora_bar_bot <- renderPlot({
    orabot <- contrast_fun()$ORA_Bot
    barplot(orabot, showCategory = 5)
  })

  output$gsea_go <- renderPlot({
    gsea_go <- contrast_fun()$GSEA_GO
    enrichplot::dotplot(gsea_go, showCategory = 5)
  })

  output$gsea_kegg <- renderPlot({
    gsea_kegg <- contrast_fun()$GSEA_KEGG
    enrichplot::dotplot(gsea_kegg, showCategory = 5)
  })
}

shinyApp(ui, server)
