## ----packages_scripts_setup,results='hide',message=F,warning=F,echo=F, error=F----
# global font setting for plots:
my.font='CMU Serif'
library(knitr)
library(parallel) # multicore processing
library(ggplot2)
library(GGally)
	theme_set(theme_bw(base_family=my.font))
	 # +
 		# 		  theme(
		#                  plot.background = element_rect(fill = "transparent", colour = NA)
 	        #           ))
 
	# a colorblind-friendly palette from
	# # http://jfly.iam.u-tokyo.ac.jp/color/
	# # substitute black with gray
	# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	# scale_color_discrete <- function(...) scale_color_manual(values=cbPalette, ...)
	# scale_fill_discrete <- function(...) scale_fill_manual(values=cbPalette, ...)
library(maps)
library(maptools)
library(mapproj)
library(qvalue)
library(pcaMethods)
library(ggbiplot)
library(scales)
library(Hmisc)
library(tidyr)
library(dplyr)
library(stringr)
# knitr defaults:
opts_chunk$set(autodep=T,
		echo=F,
		comment=NA,
		message=F,
		warning=F, 
		out.width='\\textwidth', 
		dev='cairo_pdf', # font declaration in plots only works when using cairo
		size='footnotesize') 
knit_hooks$set(crop = hook_pdfcrop)
opts_knit$set(concordance=TRUE,self.contained=TRUE) # allowing synctex between PDF and .Rnw
options(scipen = 2, digits = 2) 

## ----dom,results='hide',crop=T-------------------------------------------
load('SOM/autotyp.dom.rData')
autotyp.dom$mysymbol <- with(autotyp.dom, ifelse(sbranch %in% 'Italic', 24, ifelse(mbranch %in% 'Indo-Iranian',25,21)))
autotyp.dom$mycolor <- with(autotyp.dom, ifelse(DOM,'black','white'))
	if (packageVersion('maps') < 3) {
	map('world', projection='gall', param=0, interior=F, xlim=c(-21,177),ylim=c(7.5,73), col='darkgrey') 
	} else {
	# the new database in 'maps' makes troubles with not obeying the xlim and ylim stuff, so we use:	
	map('legacy_world', projection='gall', param=0, interior=F, xlim=c(-21,177),ylim=c(7.5,73), col='darkgrey')	
	}
with(subset(autotyp.dom, !macrocontinent %in% c('Africa','Americas') & !area  %in% 'Oceania' & latitude>7), points(mapproject(x=longitude, y=latitude), cex=.7, col='black', lwd=.6, pch=mysymbol, bg=mycolor))

## ----isolates,cache=T----------------------------------------------------
load('SOM/2015-04-01T16-01autotyp.Rdata')
autotyp.gg <- autotyp.backbone[,c(autotyp.essentials[1:16], 'genesis', 'modality')]
autotyp.gg %>% 
	filter(!genesis %in% 'creole' & !modality %in% 'signed') %>% 
	group_by(stock) %>% 
	dplyr::summarize(members=n()) %>% 
	select(members) %>% unlist -> autotyp.stock.size

## ----uniformity, cache=T-------------------------------------------------
load('SOM/typology.Rdata')
uniformity <- do.call('rbind', 
                      lapply(typology.data, function(l) {
 			     l[,1:3] %>% 
			     filter(complete.cases(.)) %>% 
			     group_by(stock) %>% 
			     summarise_each(funs(ntypes=n_distinct, length), 2)  %>% 
			     filter(length>4) %>% 
			     mutate(uniform = ntypes==1, nonuniform = ntypes > 1) %>% 
			     dplyr::summarise(
				        prop.uniform=mean(uniform), 
		             		count.uniform=sum(uniform), 
			             	prop.non.uniform=mean(nonuniform),
		             	   	count.non.uniform=sum(nonuniform),
		             	        fams=n()
					) -> x
		             data.frame(var=names(l)[2],x) %>% 
  			     filter(fams>0)
		     	     }
 			    ))
# ggplot(uniformity, aes(x=prop.uniform)) +
#  geom_histogram(aes(y=..count../sum(..count..)), fill='darkgrey') +
#  scale_y_continuous(name='Relative frequency', labels = percent_format()) +
#  xlab("Proportion of uniform families per variable\n Note: families with at least 5 members;\n 719 variables from WALS [with many recodings] and AUTOTYP") +
#  theme(axis.title.x=element_text(lineheight = 1.2, vjust=-.8))


## ----data,cache=T--------------------------------------------------------
# create area-based pseudo-groups when a family straddles the area border. See Bickel (2013, Distributional biases) for justification
autotyp.backbone$ssbranch <- as.character(autotyp.backbone$ssbranch)
autotyp.backbone$ssbranch[autotyp.backbone$language %in% c('Evenki','Solon')] <- 'Inland NW Tungusic pseudogroup'
autotyp.backbone$ssbranch[autotyp.backbone$language %in% c('Negidal','Orok')] <- 'Coastal NW Tungusic pseudogroup'
autotyp.backbone$ssbranch[autotyp.backbone$sbranch %in% 'Western Malayo-Polynesian non-clade' & is.na(autotyp.backbone$ssbranch) & autotyp.backbone$CP %in% 'CP'] <- 'Eastern WMP pseudo-group'
autotyp.backbone$ssbranch[autotyp.backbone$sbranch %in% 'Western Malayo-Polynesian non-clade' & is.na(autotyp.backbone$ssbranch) & autotyp.backbone$CP %in% 'Other'] <- 'Western WMP pseudo-group'

# fix one error in the database
autotyp.backbone$CP[autotyp.backbone$language %in% 'Bima'] <- 'CP'

# subset to N(data)>250 and N(levels)<10
typology.data.raw <- lapply(autotyp.all, function(d) {merge(d[!is.na(d[,2]),], subset(autotyp.backbone, !is.na(CP), autotyp.essentials), by='LID')})
typology.data1 <- typology.data.raw[sapply(typology.data.raw, function(l) length(unique(l$LID))>=250)]
typology.data2 <- typology.data1[sapply(typology.data1, function(l) length(unique(l[,2]))<10)]

# for S and the adnominal functions it turns out there is never a difference in how the more fine-grained (SR) and the more abstract role definitions slice up the role, so we remove them from the SR codes. 
typology.data <- typology.data2[!names(typology.data2) %in% c("autotyp.basic.sr.locus.per.language$BasicLocus.POSS", "autotyp.basic.sr.locus.per.language$SimpleLocus.POSS", "autotyp.basic.sr.locus.per.language$BasicLocus.S", "autotyp.basic.sr.locus.per.language$SimpleLocus.S")]

# Note: because we limit the data to N>=250, only nonspecial roles tend to survive in the data, e.g. U but not U-def, i.e. U is strictly the nondefinite undergoer

# finally, just for safety, recode all logical vectors as factors
typology.data <- lapply(typology.data, function(i) {i[,2] <- as.factor(i[,2]); return(i)})

# report
coverage <- range(sapply(typology.data, function(d) length(unique(d$LID))))

# save the data in two batches so we can distributed this on two machines:
save(typology.data, file='SOM/typology.Rdata')
typology.data.all <- typology.data
split <- length(typology.data)/2
typology.data <- typology.data.all[1:split]
save(typology.data, file='SOM/typology1.Rdata')
typology.data <- typology.data.all[(split+1):length(typology.data.all)]
save(typology.data, file='SOM/typology2.Rdata')

## ----tp,crop=T-----------------------------------------------------------
# only languages for which we have data:
cp.lids <- Reduce(union, mclapply(typology.data, function(d) d[,1]))
autotyp.cp <- subset(autotyp.backbone, LID %in% cp.lids)
autotyp.cp$mycolor <- with(autotyp.cp, ifelse(CP %in% 'CP','black','darkgrey'))
# autotyp.cp$mycolor <- with(autotyp.cp, ifelse(CP %in% 'CP','red','black')); alpha=.6
if (packageVersion('maps') < 3) {
	map('world', projection='cylequalarea', param=0, orientation=c(90,150,0), wrap=T, interior=F, col='darkgrey', lwd=.5)
	} else {
	map('legacy_world', projection='cylequalarea', param=0, orientation=c(90,150,0), wrap=T, interior=F, col='darkgrey', lwd=.5)
	}
with(autotyp.cp, points(mapproject(x=longitude, y=latitude), cex=.6, col=NULL, lwd=.6, pch=21, bg=mycolor))

## ----results,fig.height=2.5----------------------------------------------
tp.results1 <- read.table('SOM/results1.tab', header = T, stringsAsFactors=F)
tp.results2 <- read.table('SOM/results2.tab', header = T, stringsAsFactors=F)
tp.results.tmp <- rbind(tp.results1,tp.results2)
# note: the current results [2015-04-25] were computed when the duplicates for locus of S and POSS were still in the autotyp.all dataset, so need to take them out by hand again from the results:
tp.results <- subset(tp.results.tmp, !variable %in% c("autotyp.basic.sr.locus.per.language$BasicLocus.POSS", 		
					"autotyp.basic.sr.locus.per.language$SimpleLocus.POSS", 
					"autotyp.basic.sr.locus.per.language$BasicLocus.S", 
					"autotyp.basic.sr.locus.per.language$SimpleLocus.S")
					)
# also, the bias computation scripts were run when there was a small coding error in determining the source of some variables. This is now fixed, but we need to patch it here [2015-04-26] 
tp.results$source[grepl('per\\.lang', tp.results$variable) & tp.results$source %in% 'WALS'] <- 'AUTOTYP'

# example for presentations:
# ex <- familybias(typology.data[['autotyp.wals$DRYGEN']], r.name = 'autotyp.wals$DRYGEN' , p.names = 'CP', lapplyfunc = mclapply, B=100, family.names=c('stock','mbranch','sbranch','ssbranch','lsbranch','language'))
# ex.biases <- mean(ex) %>%
# 		filter(!majority.response %in% 'diverse')
# mosaic(xtabs(Freq~., ex.biases, drop.unused.levels=T),
# 	shade=T, gp=shading_max, indepfun = function(x) sum(x^2),
# 	split_vertical=T,
# 	varnames=F,
# 	set_labels=list(CP=c('Other','Trans-Pacific'), majority.response=c('GenN','NGen','neither dominant')),
# 	rot_labels=c(0,0),
# 	just_labels=c('center','right'),
# 	margins=c(2,1,1,8),
# 	keep_aspect_ratio=F, 
#	sub='(Dryer 2005, WALS)'
# 	)

# check correlations with sample size

# ggplot(tp.results, aes(x=n.languages,y=set.fisher.p)) + geom_point() + geom_smooth()
# ggplot(tp.results, aes(x=n.languages,y=ml.fisher.p)) + geom_point() + geom_smooth()
# ggplot(tp.results, aes(x=n.languages,y=mcmc.fisher.p)) + geom_point() + geom_smooth()

# plot results					
tp.results.long <- tp.results %>% gather(key=Method, value=p, contains('.p'))  %>% filter(is.finite(p))
levels(tp.results.long$Method) <- c('set-based','tree-based (MCMC)', 'tree-based (ML)')
ggplot(tp.results.long, aes(x=p, linetype=Method)) + 
	# geom_density() +
	geom_line(stat='density') +
	scale_x_continuous(breaks=seq(0,1,.05), labels=c('0',gsub('^0','',seq(.05,1,.05)))) +
	scale_linetype_discrete(name='Methods:') +
	theme(title=element_text(face='italic'),
		axis.title.y=element_text(vjust=1),
		axis.title.x=element_text(vjust=.6),
		legend.title=element_text(face='italic'), 
		legend.key=element_rect(color="white"), 
		legend.key.width=unit(1,"line")
		) # + theme_bw(base_size = 24) 		

## ----results-summary-----------------------------------------------------
tp.results.set <- subset(tp.results, !is.na(set.fisher.p), c('variable','set.fisher.p')) # some problems
set.q <- qvalue(tp.results.set$set.fisher.p, pi0.method="bootstrap")$qvalues
tp.results.set <- data.frame(tp.results.set,set.q)

tp.results.trees <- subset(tp.results, !is.na(mcmc.fisher.p) & !is.na(ml.fisher.p), c('variable','mcmc.fisher.p','ml.fisher.p'))  # if there are NAs (problems), they are the same for both, arising from lack of convergence in BayesTraits
mcmc.q <- qvalue(tp.results.trees$mcmc.fisher.p, pi0.method="bootstrap")$qvalues
ml.q <- qvalue(tp.results.trees$ml.fisher.p, pi0.method="bootstrap")$qvalues
tp.results.trees <- data.frame(tp.results.trees, mcmc.q, ml.q)

## ----results-details, cache=T--------------------------------------------
tp.results.sign.tmp <- subset(tp.results, set.fisher.p<.05 | mcmc.fisher.p<.05) 

# fix names and add variants
name.edits <- read.csv('SOM/variable_name_editing.csv', header=T, stringsAsFactors=F)
tp.results.sign <- merge(tp.results.sign.tmp, name.edits, by.x='variable', by.y='original')
# add the variants from WALS automatically
tp.results.sign$variant_of <- with(tp.results.sign, ifelse(variant_of %in% "" & grepl('[A-Z]{6}$',variable), gsub('(.*\\$)([A-Z]{6})$','\\2', tp.results.sign$variable), variant_of))
tp.results.sign$variant_of <- with(tp.results.sign, ifelse(variant_of %in% "" & grepl('[A-Z]{6}\\d+$',variable), gsub('(.*\\$)([A-Z]{6})\\d+$','\\2', tp.results.sign$variable), variant_of))
# add remaining variable names
tp.results.sign$new <- with(tp.results.sign, ifelse(new %in% "", gsub('(.*\\$)(.*)','\\2', tp.results.sign$variable), new))

# fix the value names that we will report
tp.results.sign$set.CP.significant.residuals <- gsub('\\n',';',tp.results.sign$set.CP.significant.residuals)
tp.results.sign$set.Other.significant.residuals <- gsub('\\n',';',tp.results.sign$set.Other.significant.residuals)
# (the following is for export only so we can edit the names in a spreadsheet)
# tp.results.sign$TransPacific <- tp.results.sign$set.CP.significant.residuals
# tp.results.sign$Other <- tp.results.sign$set.Other.significant.residuals
# tp.results.sign[,c('new','set.CP.significant.residuals', 'TransPacific','set.Other.significant.residuals','Other')]
feature.edits <- read.csv('SOM/feature_names.csv', header=T, stringsAsFactors=F)
tp.res.sign.clean <- merge(tp.results.sign, feature.edits[,c('variable','TransPacific','Other')], by='variable')
tp.res.sign.clean$set.fisher.p <- round(tp.res.sign.clean$set.fisher.p,4)
tp.res.sign.clean$mcmc.fisher.p <- round(tp.res.sign.clean$mcmc.fisher.p,4)
tp.res.sign.clean$ml.fisher.p <- round(tp.res.sign.clean$ml.fisher.p,4)

# "variant of" calculation
# if a variable has an entry but occurs only once in this table under either test threshold, we remove it:
singles <- names(table(tp.res.sign.clean$variant_of))[table(tp.res.sign.clean$variant_of)==1]
tp.res.sign.clean$variant_of[tp.res.sign.clean$variant_of %in% singles] <- ""
# number of variants set.based:
tp.res.sign.clean.set <- subset(tp.res.sign.clean, set.fisher.p<.05)
variants.set <- tp.res.sign.clean.set$variant_of[!tp.res.sign.clean.set$variant_of==""]
nr.variants.set <- sum(table(variants.set))-length(unique(variants.set))  # substract one each because one of two variants is not a variant..
# number of variants mcmc.based:
tp.res.sign.clean.mcmc <- subset(tp.res.sign.clean, mcmc.fisher.p<.05)
variants.mcmc <- tp.res.sign.clean.mcmc$variant_of[!tp.res.sign.clean.mcmc$variant_of==""]
nr.variants.mcmc <- sum(table(variants.mcmc))-length(unique(variants.mcmc))  

## ----appendix, results='asis', cache=F-----------------------------------
tp.res.sign.clean %>% 
	arrange(set.fisher.p,mcmc.fisher.p) %>% 
	mutate(source=substr(source,1,1)) %>% 
	select(new, source, n.languages, contains(".p"), TransPacific, Other, variant_of)  %>% 
	latex(., file="",
	caption='\\textbf{Appendix:} Variables with a significant difference in family biases, with \\emph{p} < .05 for at least the set-based or the MCMC tree-based estimates (ordered by the the set-based \\emph{p}-values). Plus and minus signs indicate Pearson residuals that are significant under a permutation test with a .05 rejection level. Missing \\emph{p}-values mean that family bias estimates failed for statistical or computational reasons. A = \\textsc{autotyp}, W = \\textsc{wals}',
	colheads=c('\\textit{Variable}','\\textit{Source}','\\textit{N}(lang.)','\\textit{p}(sets)','\\textit{p}(MCMC)', '\\textit{p}(ML)', '\\textit{Trans-Pacific}', '\\textit{Other}', '\\textit{Variant of}'),
	col.just=c("l","l", "r", "r", "r", "r", "l", "l", "l"),
	collabel.just=c("l","l", "r", "r", "r", "r", "l", "l", "l"),
	rowname=NULL,
	booktabs=T,
	where='p',
	size='scriptsize',
	landscape=T, longtable=T, lines.page=4000, # lines.pages= needed because the default inserts too many breaks
	label='tab-details'
	) -> x

