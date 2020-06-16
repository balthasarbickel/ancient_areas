# Family Bias test of the Trans-Pacific area

# This file has '1' at the end because it works on the first half of the data, stored in typology1.rData

# Load the nessesary packages
library(devtools)
library(ape)
library(geiger)
library(parallel)
library(vcd)
source_url('https://github.com/IVS-UZH/phylo-convert/raw/master/phylo-convert.R') 
source_url('https://github.com/IVS-UZH/familybias/raw/master/R/familybias.R') 

# NOTE on run time: if the data is divided in two equal sets, each set takes about 2h30' on a 2013 6-core MacPro and nearly 3h on a 12-core 2011 PowerPC

buildResults <- function(current.df, varname, fbias.pre.correction, fbias.ml.correction, fbias.mcmc.correction)
{
	# prepare the template
	results <- data.frame(
			variable=varname,
			raw.sample.size=nrow(current.df),
			n.languages=length(unique(current.df$LID)),
			n.levels=length(unique(current.df[[varname]])),
			source=if(grepl('wals\\$',varname)) 'WALS' else 'AUTOTYP',
			set.fisher.p=NA,
			set.CP.significant.residuals=NA,
			set.Other.significant.residuals=NA,
			mcmc.fisher.p=NA,
			mcmc.CP.significant.residuals=NA,
			mcmc.Other.significant.residuals=NA,
			ml.fisher.p=NA,
			ml.CP.significant.residuals=NA,
			ml.Other.significant.residuals=NA,
			problems=NA,
			stringsAsFactors=F)
	
	compute_stats_from_fbias <- function(fbias)
	{
		fbias.m <- subset(mean(fbias), !majority.response %in% "diverse")
		fbias.m$Freq <- round(fbias.m$Freq)
		
		
		bias.df <- as.data.frame(xtabs(Freq ~ CP + majority.response, fbias.m, drop.unused.levels=T))
		bias.df <- as.matrix(xtabs(Freq~., bias.df,drop.unused.levels=T))
		
		bias.df <- bias.df[, !apply(bias.df, 2, function(x) all(x==0)), drop=T]			
		
		ctest <- coindep_test(bias.df)
				
		critical.pearson <- ctest$qdist(0.95)
		
		if(STAT_DEBUG)
		{
			cat('================\n')
			print(bias.df)
			print(ctest$expected)
			print(ctest$residuals)
			cat('Critical residual level at', critical.pearson, '\n')
		}
		
		
		get_significant_residuals <- function(residuals)
		{
			if(is.na(critical.pearson)) return(NA)
			
			x <- character(length(residuals))
			x <- ifelse(residuals >= critical.pearson, paste0('+', colnames(bias.df)), x)
			x <- ifelse(residuals <= -critical.pearson, paste0('-', colnames(bias.df)), x)
			x <- x[!x=='']
			paste0(x, collapse='\r')
		}
		
		list(
			fisher.p = fisher.test(bias.df)$p.value, 
			CP.significant.residuals = get_significant_residuals(ctest$residuals['CP', ]),
			Other.significant.residuals = get_significant_residuals(ctest$residuals['Other', ])
		)
	}
		
	# add the results
	if(!is.null(fbias.pre.correction))
	{		
		stats <- compute_stats_from_fbias(fbias.pre.correction)
		
		results$set.fisher.p <- stats$fisher.p
		results$set.CP.significant.residuals <- stats$CP.significant.residuals
		results$set.Other.significant.residuals <- stats$Other.significant.residuals
	}

	if(!is.null(fbias.ml.correction))
	{
		stats <- compute_stats_from_fbias(fbias.ml.correction)		
		
		results$ml.fisher.p <- stats$fisher.p
		results$ml.CP.significant.residuals <- stats$CP.significant.residuals
		results$ml.Other.significant.residuals <- stats$Other.significant.residuals
	}
	
	if(!is.null(fbias.mcmc.correction))
	{		
		stats <- compute_stats_from_fbias(fbias.mcmc.correction)
		
		results$mcmc.fisher.p <- stats$fisher.p
		results$mcmc.CP.significant.residuals <- stats$CP.significant.residuals
		results$mcmc.Other.significant.residuals <- stats$Other.significant.residuals
	}
	
	
	
	
	results
}


# This is the outer interface to bias estimation
# Does the preprocessing and generates the report

processVariable <- function(current.df, varname)
{
	# we need to ensure that there are no duplicate names in the taxonomy. Best do
	# it now, so that things can be consistent. I will just go radical
	# here and rename ALL taxa, this is the simplest way. Also, I will make sure
	# that all taxa names are compatible with whatever BayesTraits wants. 
	
	for(l in taxonomyLevels[-length(taxonomyLevels)])
	{
		current.df[[l]] <- ifelse(!is.na(current.df[[l]]), paste(current.df[[l]], l), NA)
	}
	
	# check the tips
	lowest_taxon <- taxonomyLevels[[length(taxonomyLevels)]]
	
	# make sure that there are no NAs on the lowest level
  stopifnot(!any(is.na(current.df[[lowest_taxon]])))

	# similarly, make sure that every tip has a unique name and that there are no
	# spaces in the language names as this appears to throw BayesTraits off	
	current.df[[lowest_taxon]] <- gsub(' +', '', current.df[[lowest_taxon]])
	
	# check if there are duplicates on the lowest level, if yes, we want to
	# introduce subsystems
	if(any(duplicated(current.df[[lowest_taxon]])))
	{
		lowest_taxon1 <- paste0(lowest_taxon, '-subsystem')
		current.df[[lowest_taxon1]] <- factor(paste(current.df[[lowest_taxon]], 1:nrow(current.df), sep='-'))
		attr(current.df, 'taxa') <- c(taxonomyLevels, lowest_taxon1)
	}
	else attr(current.df, 'taxa') <- taxonomyLevels
			
	# the inner processor might crash
	# if it crashes, we simply report an error		
	r <- buildResults(current.df, varname, NULL, NULL, NULL)
	
	result <- tryCatch(
		{
		  # run the family bias estimation (before BayesTraits correction)
	    fbias.pre.correction <- familybias(current.df, family.names = attr(current.df, 'taxa'), r.name = varname, p.names = 'CP', lapplyfunc = mclapply, B=B)
		
		  r <- buildResults(current.df, varname, fbias.pre.correction, NULL, NULL)
		  .processVariable(current.df, varname, fbias.pre.correction)
		},
		error = function(e)
		{
			r$problems <- as.character(e$message)
			
			r
		})
		
	result	
}


# The internal workhorse
.processVariable <- function(current.df, varname, fbias.pre.correction)
{
 # prepare for BayesTraits correction
 bayes_traits_data <- current.df
 names(bayes_traits_data)[match(varname, names(bayes_traits_data))] <- '.VARIABLE.'
 # prepare the list of units we can submit to BayesTraits, i.e. units with non-unform value distributions:
 units <- droplevels(subset(fbias.pre.correction$large.families.estimate, majority.prop<1 | is.na(majority.prop), select=c('family.name','taxonomic.level')))

 # run BayesTraits and gather the stats
 bayes_traits_output <- .runBayesTraitsAndParseResults(bayes_traits_data, units)
 
 # build the override tables for the ML model
 bias_overrides_ml <- lapply(names(bayes_traits_output), function(unit.name) {
	 # extract the stats
	 stats <- bayes_traits_output[[unit.name]]
	 # the family is biased of if the ML p value is < 0.1
	 # then we take the majority values computed by BayesTraits
	 if (stats$lr.p.ml < 0.1) 
			stats$maj.value.data.ml
	 else
			list(mj.prop=NA, mj.val=NA)
 })
 names(bias_overrides_ml) <- names(bayes_traits_output)

 # build the override tables for the MCMC model
 bias_overrides_mcmc <- lapply(names(bayes_traits_output), function(unit.name) {
	 # extract the stats
	 stats <- bayes_traits_output[[unit.name]]
	 # the family is biased of if the bayes factor is >2
	 # then we take the majority values computed by BayesTraits
	 if ((stats$lr.mcmc > 2)) 
		 stats$maj.value.data.mcmc
	 else
		 list(mj.prop=NA, mj.val=NA)
 })
 names(bias_overrides_mcmc) <- names(bayes_traits_output)
 
 # run the family bias estimation (with BayesTraits ML correction)
 fbias.ml.correction <- familybias(current.df, family.names = attr(current.df, 'taxa'), r.name = varname, p.names = 'CP', lapplyfunc = mclapply, B=B, bias.override = bias_overrides_ml)

 # run the family bias estimation (with BayesTraits MCMC correction)
 fbias.mcmc.correction <- familybias(current.df, family.names = attr(current.df, 'taxa'), r.name = varname, p.names = 'CP', lapplyfunc = mclapply, B=B, bias.override = bias_overrides_mcmc)
 
 buildResults(current.df, varname, fbias.pre.correction, fbias.ml.correction, fbias.mcmc.correction)	
}



# runs the BayesTraits algorithm for the data
.runBayesTraitsAndParseResults <- function(bayes_traits_data, units)
{	
	# prepare the trees
	suppressWarnings({
		trees <- as.phylo(bayes_traits_data, levels=attr(bayes_traits_data, 'taxa'), roots=units)
		if(!'multiPhylo' %in% class(trees)) trees <- list(trees)
	})
	
	# init the path
	unlink(bayesTraitsDataPath,recursive=T)
	dir.create(paste(bayesTraitsDataPath, 'results', sep='/'), recursive=T, showWarnings=F)
	
	# write out the tree and data for every unit
	# run BayesTraits and parse the input back
	mclapply(seq_along(trees), function(i)
	{
		# extract the appropriate data 
		current.tree <- trees[[i]]
		unitname <- current.tree$node.label[1]
		lowest_taxon <- attr(bayes_traits_data, 'taxa')[[length(attr(bayes_traits_data, 'taxa'))]]
		current.data <- bayes_traits_data[bayes_traits_data[[lowest_taxon]] %in% current.tree$tip.label, c(lowest_taxon, '.VARIABLE.')]

		# format the data so that BayesTraits can read it
		# there can be no node labels in trees
		current.tree$node.label <- NULL
		# and the trees should be cleaned up
		current.tree <- collapse.singles(current.tree) 
		# recode the variable levels as letters
		VARIABLE.recoding <- data.frame(original = unique(current.data$.VARIABLE.))
		VARIABLE.recoding$recoded <- factor(LETTERS[1:nrow(VARIABLE.recoding)])
	  current.data$.VARIABLE. <- VARIABLE.recoding$recoded[match(current.data$.VARIABLE., VARIABLE.recoding$original)]

    # output the files
    output_file_path <- paste0(bayesTraitsDataPath,'/','UNIT', i)
		write.nexus(current.tree, file=paste0(output_file_path, ".tree"))
		write.table(current.data, file=paste0(output_file_path, ".data"), row.names = F, quote = F, col.names = F, sep = "\t")
		
		# write out the command files
		cfile <- function(.) paste0(bayesTraitsDataPath, '/commands.UNIT', i, '.', ., '.txt')
		write(c('1','1','RestrictAll qAB','run'),sep='\n', file=cfile('ER'))
		write(c('1','2','RestrictAll qAB','run'),sep='\n', file=cfile('ER.MCMC'))
		write(c('1','1','run'),sep='\n', file=cfile('ARD'))
		write(c('1','2','run'),sep='\n', file=cfile('ARD.MCMC'))	
   	 
		# run all the commands and parse the output		
		rfile <- function(.) paste0(bayesTraitsDataPath, '/results/result.UNIT', i, '.', ., '.txt')
		output <- list()
		for(command in c('ER', 'ER.MCMC', 'ARD', 'ARD.MCMC'))
		{	
			exe <- paste('./BayesTraitsV2', 
			             paste0(output_file_path, ".tree"),
									 paste0(output_file_path, ".data"), 
									  '<', cfile(command), '>', rfile(command))
			system(exe)
			
			tryCatch(
				output[[command]] <- parseBayesTraitOutput(rfile(command)),
				error = function(e)
				{
					stop('Unable to parse BayesTraits results')
				})
			
	  }		
		
		# helper func to determine the majority value for the BayesTraits results
		find_majority_value <- function(tab) 
		{
			# get the transition probabilities
			# its the last row of the table and the columns are on form qXY
			trans_prob <- unlist(tab[nrow(tab), grepl('^q[A-Z][A-Z]$', names(tab)), drop=T])			
			# we want to look for the highest transition rate
			max.transition.rate <- names(trans_prob)[which.max(trans_prob)]
			
			# and the majority value is last letter
			maj.val <- substr(max.transition.rate, nchar(max.transition.rate), nchar(max.transition.rate))
			
			# recode it back and return it together with the proportion
			maj.val.original <- as.character(VARIABLE.recoding$original[match(maj.val, VARIABLE.recoding$recoded)])			
			list(mj.prop=sum(as.character(current.data$.VARIABLE.) %in% as.character(maj.val))/nrow(current.data), mj.val=maj.val.original)
		}
		
		
		# gather the result statistics
		stats <- list(
			# 1. the ML likelihood ratio and majority value of the ARD model
			lr.ml=2*(output$ARD$Lh-output$ER$Lh),
			maj.value.data.ml = find_majority_value(output$ARD),
			# 2. the bayes factor of MCMC evaluation and the majority value of the 
			# ARD.MCMC model
			lr.mcmc=2*(output$ARD.MCMC$Harmonic.Mean[nrow(output$ARD.MCMC)]-output$ER.MCMC$Harmonic.Mean[nrow(output$ER.MCMC)]),
			maj.value.data.mcmc = find_majority_value(output$ARD.MCMC)
		)	
		stats$lr.p.ml <- 1 - pchisq(stats$lr.ml, 2*choose(length(unique(current.data$.VARIABLE.)),2)-1)
		
		
		# return the stats
		stats 
	}) -> results
	
	names(results) <- sapply(trees, function(.) .$node.label[1])
	
	# check if there are any errors
	errors <- sapply(results, inherits, 'try-error')
	if(any(errors))
	{
		errors <- results[errors]
		errors <- sapply(errors, function(.) attr(., 'condition')$message)
		errors <- paste(names(errors), errors, sep=':', collapse='; ')
		stop(errors)
	}
	
	
	results
}


parseBayesTraitOutput <- function(path)
{
	lines <- readLines(path)
	
	# extract only the table 
	# we use a very simple heuristic here: tables contain at least 5 tabs
	tbl <- lines[grepl('(\t[^\t]*){5,}', lines, perl=T)]
	tbl <- gsub('[[:space:]]+$', '', tbl)
	
	read.table(textConnection(tbl), header=T, sep='\t')
}


# Load the data sets
load('typology1.rData')


# Taxononomy levels
taxonomyLevels = c("stock", "mbranch", "sbranch", "ssbranch", "lsbranch", "language")

# how many cores?
options(mc.cores = 10)

# how many interpolations
B = 1000

bayesTraitsDataPath = 'bayes_traits_data'


# the results go into this file
outputPath = 'results1.tab'

# st this to true to debug the statistics
STAT_DEBUG=FALSE



for(i in 1:length(typology.data))
#for(i in 7)
{
	# get the variable name
	varname <- names(typology.data)[i]
	
	# process the variable
	timing <- system.time({
	 	result <- processVariable(unique(typology.data[[varname]]), varname)
	})
	
	cat(paste(varname, 'processed in', timing[3], 'sec.'))
	if(!is.na(result$problems)) cat(paste(' Problems:', result$problems))
	cat('\n\n')
		
	if(i == 1)
		write.table(result, file = outputPath, row.names=F, col.names=T, sep='\t')
	else
		write.table(result, file = outputPath, row.names=F, col.names=F, append=T,sep='\t')
}