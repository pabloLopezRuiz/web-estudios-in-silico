# Proyecto de creación de una web para la visualización de resultados del experimento asociado a GSE124799

En este proyecto he intentado reproducir un análisis de expresión diferencial
sobre los datos de GSE124799. Donde se comparaba el transcriptoma de células
de raton en 3 condiciones distintas: ejercico, sedentarismo y sedentarismo
tras unas semanas de ejercico. Para hacer una visualización interactiva de
los datos he creado una web app. Lamentablemente, tras intentar subir la web
al repositorio de shiny, uno de los objetos que necesito cargar es muy pesado
y por tanto la app se bloquea, por esto adjunto la app como fichero de R.

# Estructura del directorio

#### **web_esb/app.R**
Archivo de la app. Para visualizarla se debe de ejecutar este archivo.


#### **data**
Fichero con archivos generados en el análisis.

#### **GSE124799**
Fichero con el SummarizedExperiment archivo (*gse.RData*), además contiene
los archivos generados en el contraste.

#### **processing.R**
Archivo donde se transforma de los datos crudos subidos a GEO y los metadatos
al archivo Summarized Experiment usado en el resto del análisis.

#### **DEA.R**
Archivo del análisis de expresión diferencial, con algunos gráficos generados

#### **ORA.R**
Archivo del enriquecimiento funcional, contiene un ORA, un GSEA y un análisis
de STRING.

#### **pipeline_completo.qmd\html**
Archivo con un documento generado con todo el análisis hecho en processing.R, 
DEA.R y ORA.R comprimido en un único informe con los pasos detallados.


