shinyServer(function(input, output, session) {
	obj <- reactiveValues() # we work with this data!

	# Data #
	obj$importedData <- NULL
	obj$importedFeatureSummary <- NULL
	obj$subsetData <- NULL
	obj$subsetFeatureSummary <- NULL
	obj$transformedData <- NULL
	obj$transformedFeatureSummary <- NULL
	obj$is_it_gene_expression_data <- NULL
	obj$seuratObj <- NULL

	obj$transform_method <- NULL
	# obj$clearDataFile <- TRUE
	obj$file_data_csv <- NULL
	obj$file_data_rdata <- NULL
	obj$press_data_tabPanel <- "Display data"
	obj$imported_data_or_file_name <- NULL

	

	# Causality Test #
	obj$allResultList <- NULL

	obj$kegg_gene_summary <- NULL
	obj$runCausalTestQ <- FALSE
	obj$clearDAGFile <- TRUE
	obj$amatCandidateList <- NULL
	obj$candidateDAG_df <- NULL
	obj$amatTestedIndex <- NULL
	obj$amatTestedList <- NULL
	obj$amatTestedIndex_valid <- NULL
    obj$amatTestedList_valid <- NULL
	obj$networkVisualization_input_DAG <- NULL
	obj$DAG_test_currnet_visnetwork <- NULL
	obj$download_DAG_test_flag <- NULL

	obj$dag_test_input_mode <- NULL
	obj$CIT_method <- NULL
	obj$alpha <- NULL
	obj$num_of_tested_dag <- NULL
	obj$causal_test_datasets <- NULL


	observeEvent(input$btn_home_click_me, {
		updateTabsetPanel(
			session = session, 
			inputId = "navbar", 
			selected = "User manual in Chinese"
		)
	})
	
	source("./server/data.R", local = TRUE)
	source("./server/causality.R", local = TRUE)
	# source("./server/gene_information.R", local = TRUE)
	source("./server/r_sessionInfo.R", local = TRUE)
})