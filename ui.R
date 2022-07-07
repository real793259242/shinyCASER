shinyUI(
    fluidPage(
        # rclipboard::rclipboardSetup(),
        includeCSS("www/styles.css"), 
        shinyjs::useShinyjs(), 
        # theme = shinythemes::shinytheme("flatly"), 

        navbarPage(
            title = "shinyCASER",
            id = "navbar", 
            windowTitle = "shinyCASER",
            #fluid = FALSE,
            # theme = shinythemes::shinytheme("flatly"), 
            # footer = tagList(
            #     p("Author: Jielong Huang")
            # ),
            # shinyjs::useShinyjs(), 
            # tags$head(
            #     tags$link(rel = "stylesheet", 
            #               type = "text/css", 
            #               href = "styles.css")
            # ), 

            #####################
            # Home tab
            #####################
            tabPanel(title = "Home", 
                    icon = icon("home"),
                    # shinyLP::jumbotron(
                    #     "Welcome to shinyCASER!", 
                    #     p(HTML("This is an online <b>shiny</b>-based application for exploring <b>F</b>eature <b>A</b>ssociation and <b>C</b>ausali<b>T</b>y.")),
                    #     buttonLabel = "Click Me")
                    div(class = "jumbotron", 
                        h1("Welcome to shinyCASER!"), 
                        tags$p(HTML("This is an online <b>shiny</b>-based application for <b>CA</b>u<b>S</b>al <b>E</b>n<b>R</b>ichment analysis.")), 
                        tags$p(HTML("Click <b>Manual</b> button to learn how to use shinyCASER.")), 
                        myActionButton(
                            inputId="btn_home_click_me", 
                            label="Manual", 
                            btn.style = "btn-primary btn-lg"
                        )
                    ), 
                    fluidRow(
                        column(6, 
                            shinyWidgets::panel(
                                heading = "Features", 
                                tags$p("")
                            )
                        ), 
                        column(6, 
                            tags$h3(icon("id-card"), "Application developers"), 
                            tags$a("Jielong Huang", href = "https://github.com/real793259242"),
                            tags$br(), 
                            tags$h3(icon("envelope"), "Contact us"), 
                            tags$p(HTML("Any bugs or suggestions please feel free to send to: <a href = \"mailto: zhuhao@smu.edu.cn\">zhuhao@smu.edu.cn</a>."))
                        )
                    )
            ),  

            #####################
            # Data tab
            #####################
            tabPanel(
                title = "Data", 
                icon = icon("table"), 
                # value = "table",
                uiOutput("ui_data")
            ),  


            ###############
            # Causality tab
            ###############
            tabPanel(
                title = "Step 1 - Input DAGs", 
                uiOutput("ui_causal_test_step1")
            ), 


            tabPanel(
                title = "Step 2 - Run causality test", 
                uiOutput("ui_causal_test_step2")
            ), 


            ##########
            # Help tab
            ##########
            navbarMenu(
                title = "", 
                icon = icon("question-circle"),
                tabPanel("User manual", 
                       icon = shiny::icon("book")
                ), 
                tabPanel("User manual in Chinese", 
                       icon = shiny::icon("book"), 
                       fluidPage(
                           titlePanel("User manual in Chinese"), 
                           fluidRow(
                               column(12, htmltools::tags$iframe(src = "manual_Chinese.html", width = '100%', height = 1000, style = "border:none;"))
                           )
                       )
                ), 
                tabPanel("Q & A", 
                       icon = shiny::icon("question"), 
                       fluidPage(
                           titlePanel("Q & A"), 
                           # https://community.rstudio.com/t/including-a-html-file-rendered-from-rmarkdown-in-r-shiny-apps-using-shinydashboard-is-causing-tabitems-to-break/84819
                           fluidRow(
                               column(12, htmltools::tags$iframe(src = "QandA.html", width = '100%', height = 1000, style = "border:none;"))
                           )
                       )
                ), 
                tabPanel("R session information",
                       icon = shiny::icon("r-project"), 
                       fluidPage(
                           titlePanel("R session information"), 
                           verbatimTextOutput("sessionInfo_text") %>% withSpinner(type = 6)
                       )
                )
            )
        )
    )
)