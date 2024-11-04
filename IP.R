# Load or installed the packages required
required_packages<- list('dplyr',
	'BiocManager',
	'Biostrings',
	'doParallel',
	'openxlsx',
	'stringr',
	'msa')
lapply(required_packages, function (pack)
{
	if (!require(pack, quietly=T, character.only = T))
	{
		install.packages(pack, repos = "https://cloud.r-project.org")
	}
})
lapply(required_packages, function (pack)
{
	if (!require(pack, quietly=T, character.only = T))
	{
		BiocManager::install(pack)
	}
})

# IsoPrimer must run within its folder
# or temp file removal may cause damage
if (!grepl('IsoPrimer', getwd()))
{
	stop('IsoPrimer must run from within its folder', call.=FALSE)
}

# Read in the command line arguments
args<- commandArgs(trailingOnly=TRUE)

# Retrieve the command line arguments and set them as variables
# for the different pipeline steps

if (length(args)!=10)
{
	m0<-"The following arguments must be supplied:
		\t\t<threads number>
		\t\t<kallisto>
		\t\t<primer3>
		\t\t<primersearch>
		\t\t<expression threshold %>
		\t\t<primersearch mismatch %>
		\t\t<transcriptome>
		\t\t<annotation>
		\t\t<genome>
		\t\t<debug>"
	system(paste('echo [IsoPrimer] $(date) "', as.character(m0), '" > IP_Log.out'))
	stop(call.=FALSE)
}

THREADS <- args[1]
KALLISTO <- as.character(args[2])
PRIMER3 <- as.character(args[3])
PRIMERSEARCH <- as.character(args[4])
PEXPRESSION <- as.character(args[5])
PMISMATCH <- as.character(args[6])
TRANSCRIPTOME <- as.character(args[7])
ANNOTATION <- as.character(args[8])
GENOME <- as.character(args[9])
DEBUG <- !eval(parse(text=args[10]))
Sys.setenv(THREADS = THREADS)
Sys.setenv(KALLISTO = KALLISTO)
Sys.setenv(PRIMER3 = PRIMER3)
Sys.setenv(PRIMERSEARCH = PRIMERSEARCH)
Sys.setenv(PMISMATCH = PMISMATCH)
Sys.setenv(TRANSCRIPTOME = TRANSCRIPTOME)
Sys.setenv(ANNOTATION = ANNOTATION)
Sys.setenv(GENOME = GENOME)
m1<-paste0("Arguments supplied: \n",
	"\t\tKallisto path: ", KALLISTO, "\n",
	"\t\tPrimer3 path: ", PRIMER3, "\n",
	"\t\tPrimerSearch path: ", PRIMERSEARCH, "\n",
	"\t\tExpression threshold: ", paste0(PEXPRESSION, '%'), "\n",
	"\t\tPrimerSearch mismatch: ", paste0(PMISMATCH, '%'), "\n",
	"\t\tThreads number: ", THREADS, "\n",
	"\t\tTranscriptome (.fa): ", TRANSCRIPTOME, "\n",
	"\t\tAnnotation (.gtf): ", ANNOTATION, "\n",
	"\t\tGenome (.fa): ", GENOME, "\n")
system(paste0('echo [IsoPrimer] $(date) Primer design started: "', as.character(m1), '"', ' > IP_Log.out'))

# Transcript quantification
# Build the Kallisto index and start the quantification process

if (!file.exists('mod_transcriptome.fa'))
{
	system('bash headmod.sh || echo "[IsoPrimer.header_modder.sh] $(date) Please modify the headers of the transcriptome FASTA file" >> IP_Log.out')
	if (!file.exists('mod_transcriptome.fa')) { stop('headmod.sh could not create mod_transcriptome.fa, please modify the headers of the transcriptome FASTA file', call.=FALSE) }
}
TRANSCRIPTOME <- paste0(getwd(), '/', 'mod_transcriptome.fa')
Sys.setenv(TRANSCRIPTOME = TRANSCRIPTOME)
system('cd quantification/ && bash kallister.sh && echo [IsoPrimer.Kallisto] $(date) Running Kallisto >> ../IP_Log.out')
system('cd quantification/KA_CountingOutput/ && bash kalcounting.sh && echo [IsoPrimer.Kallisto] $(date) Quantification complete >> ../../IP_Log.out')
kal<-read.table('quantification/KA_CountingOutput/kalcounts.tsv', header=T, sep='\t', as.is=T)
if (length(which(duplicated(kal$t_name)))!=0) {kal<-kal[-which(duplicated(kal$t_name)),]}
row.names(kal)<-kal$t_name
kal<-kal[, -1, drop=F]

# cln stores info about the driver lncRNA isoforms
# tplus also, but sjs positions are numbered and not identified by colons
# everything is in wsp.RData which is the final product of prepper.R

system('echo [IsoPrimer.prepper.R] $(date) Annotation parsing started >> IP_Log.out')
ifelse(file.exists('wsp.RData'), load('wsp.RData'), source('prepper.R'))
system('echo [IsoPrimer.prepper.R] $(date) Annotation parsing complete >> IP_Log.out')

vercutter <- function (ids)
{
	return(sapply(strsplit(as.character(ids), '\\.'), `[`, 1))
}

tplus[which(tplus$exons_junction==-2),4]<-''
length(unique(tplus[,1]))

# let's choose lncRNAs: the isoform with the highest TSL AND 
# "common sjs"
# each sj has its position there's no such a thing as common sj
# score sjs according to sj-5:sj+5 snippets that are common to most isoforms
swag <- do.call(rbind, lapply(unique(tplus$name), function (e) 
{
	thr <- 5
	# For testing purposes change e in
	# ENSG00000224609 for FGGY ENSG00000227811 for INKA2-AS1 and ENSG00000225937 for PCA3
	hok <- which(cln[,1]==e)
	system(paste('echo [IsoPrimer] $(date) now processing', e, '>> IP_Log.out'))
	pol <- cln[hok,5]
	# 'pul' is a list of numeric vectors: each ordered vector reports the
	# occurrence of a sj across all isos.
	#
	# Primer3 is allowed to draw primers on every sj
	# if genes have only 1 transcript
	# in this case, 'pol' becomes a vector of increasing integers
	# as long as the number of sjs in the isoform
	#
	# if there are no sjs on that one single isoform,
	# it draws primers anywhere 
	# a non-existent chosen sj gets invoked
	# and an empty field is passed to Primer3
	# the program draws everywhere
	if (length(pol)==1)
	{
		hau <- length(which(strsplit(pol[1], "")[[1]]==':'))
		ifelse(hau!=0, pul <- list(seq(1, hau)), pul <- list(1))
	} else
	{
		pul <- lapply(pol, function (j) 
		{
			dec <- strsplit(j, "")[[1]]
			idx <- which(dec==':')
			if (length(idx)==0) {ex <- 0} else
			{
				ex <- as.vector(sapply(idx, function (i)
				{
					if (i<=thr) {thr <- i-1}
					if (is.na(i+thr)) {thr <- length(dec)-i}
					word <- paste0(dec[(i-thr):(i+thr)], collapse='')
					length(grep(word, pol, fixed=T))
				}))
			}
		})
	}
	levs <- as.integer(levels(factor(unlist(pul))))
	levs <- levs[order(levs, decreasing=T)]
	# Make core clusters according to the THREADS option
	# and pass the sjs (levs) in chunks
	cores <- as.integer(THREADS)
	if (cores<length(levs))
	{
		chunks <- length(levs)%/%cores
	} else
	{
		cores<- length(levs)
		chunks<- 1
	}
	cl <- makeCluster(cores)
	registerDoParallel(cl, cores=cores)
	swagger <- do.call(rbind, lapply(seq(chunks), function (chunk)
	{
		# The number of cores is the same as the length
		# the levs vector is split in
		ifelse(chunks!=1, cope<- cores-length(split(levs, ceiling(seq_along(levs)/(chunks+1)))[[chunk]]), cope<-0)
		swaggerino <- foreach(w=((chunk-1)*cores+1):((chunk*cores)-cope),
			.combine='rbind',
			.export=c('levs', 'pul', 'hok', 'kal', 'cln', 'tplus', 'transcripts', 'TRANSCRIPTOME', 'PRIMER3', 'PRIMERSEARCH', 'PEXPRESSION'),
			.packages='dplyr') %dopar%
		{
			Sys.setenv(LVL = w)
			Sys.setenv(TRANSCRIPTOME = TRANSCRIPTOME)
			Sys.setenv(PRIMER3 = PRIMER3)
			Sys.setenv(PRIMERSEARCH = PRIMERSEARCH)
			mx <- levs[w]
			chosed <- lapply(pul, function (v) {which(v==mx)})
			chosen <- unlist(lapply(chosed, function (v)
				{
					if (length(v)>0) {out <- paste0(v, collapse=' ')} else {out <- 0}
				}), use.names=F)
			(yolo <- as.data.frame(cbind(cln[hok, 3], chosen)))
			names(yolo)<-c('t_name', 'chosen_sj')
			tminus <- right_join(tplus, yolo, by='t_name')
			# Add transcript support level
			# Choose transcripts with TSL≤3 OR NA with the chosen sj
			tminus<-left_join(tminus, transcripts[,2:3], by='t_name')
			tminus[,5:6]<-apply(tminus[,5:6], 2, as.character)
			tminus[which(tminus$t_support_lvl=='NA'),6]<-0
			tminus<-tminus[which(tminus$t_support_lvl<=3&tminus$chosen_sj!="0"),]
			# only isos with the chosen sj will be passed to primer3
			plc<-lapply(strsplit(as.character(tminus$exons_junction), ' '), as.numeric)
			hld<-lapply(strsplit(as.character(tminus$chosen_sj), ' '), as.numeric)
			tminus['chosen_sj']<-unlist(mapply(function(v, k) {paste0(v[c(k)], collapse=' ')}, plc, hld), use.names=F)
			# the following line introduces tollerance towards
			# isoforms without sj (compare above)
			tminus[which(tminus$chosen_sj=='NA'),5]<-''
			tminus<-tminus[,c(1:3, 5)]
			if (nrow(tminus)==0) { spl <- setNames(data.frame(matrix(ncol=4, nrow=0)), c('name', 't_name', 'forward', 'reverse'))} else
			{
				# save a filtered sequence/sj list
				# then draw primers with primer3
				write.table(tminus,
					paste0(w, '_transcript_seqs_and_junctions.txt'),
					col.names=F, row.names=F, sep='\t', quote=F, na='')
				system('bash ppicker.sh', ignore.stdout=T)
				system('sleep 20')
				# Cleanup the primer3 output
				spl <- read.table(paste0(w, '_primers.txt'),
					 as.is=T, header=F, sep='\t')
				spl <- spl[,1:4]
				names(spl) <- c('name', 't_name', 'forward', 'reverse')
				spl <- spl[which(!is.na(spl$forward)),]
				# Remove duplicate couples
				# to cut primersearch time
				spl <- spl[!(duplicated(paste(spl$forward, spl$reverse))),]
			}
			# Check if any primer passes filters
			# return a blank table otherwise
			if (nrow(spl)==0) { spl[,] } else
			{
				# Save the primer list and run primersearch
				write.table(spl, paste0(w, '_targets.txt'), col.names=F, row.names=F, sep='\t', quote=F)
				system('bash psearcher.sh', ignore.stdout=T)
				system('sleep 20')
				# Parse the primersearch output and get statistics
				system('bash psearch_reshaper.sh')
				system('sleep 20')
				# Read the output of primersearch
				# and ∀ lncRNA, grade primers
				psc <- read.table(paste0(w, '_psrch_resh.txt'), header=F, sep='\t', as.is=T)
				names(psc) <- c('symbol', 'coverage', 'too_long', 'detail')
				psc[which(is.na(psc$too_long)),3]<-''
				muh <- unique(spl$name)
				tot <- length(transcripts[which(transcripts$name==muh),2])
				# Some amplification products ≥ max specified length in rock.bakup +25 bp?
				# in case, couples are discarded
				too <- psc$too_long==''
				# What's the maximum difference between the longest and shortest amplicons?
				# A penalty is calculated for couples that
				# generate different amplicons 
				# (the penalty is proportional to the length difference).
				# Primers making amplicons closer to the average of the lenght interval specified by the user
				# are preferred
				charlie <- do.call(rbind, lapply(psc$detail, function (q)
				{
					bulba <- strsplit(q, '[ _]')[[1]]
					saur <- which((seq_len(length(bulba)) %% 2)==0)
					pika <- as.integer(bulba[saur])
					chu <- which((seq_len(length(bulba)) %% 2)==1)
					#ash <- mean(eval(parse(text= system('cat rock.bakup | grep PRIMER_PRODUCT_SIZE_RANGE | cut -d= -f2 | tr "-" ":"', intern=T))))
					#c(max(pika)-min(pika), paste(bulba[chu], collapse=' '), (abs((ash-min(pika))/1000)))
					c(max(pika)-min(pika), paste(bulba[chu], collapse=' '), (abs((200-min(pika))/1000)))
				}))
				# the maximum acceptable length difference 
				# between amplicons is 25 nucleotides
				golf <- as.integer(charlie[,1])<25
				# Discard a couple if there's an aspecific amplification product
				# or one of the previous criteria was not respected
				par <- sapply(1:length(psc$coverage), function (j)
				{
					chop <- strsplit(psc$coverage[j], ' ')[[1]]
					chk <- chop==muh
					iso <- length(which(chk==T))
					asp <- length(which(chk==F))
					asp <- asp==0
					tool <- too[j]
					delta <- golf[j]
					if (asp & tool & delta) {print(iso)} else {print('not_eligible')}
				})
				# favor a primer couple
				# if it covers an expressed isoform
				# according to KALLISTO's estimations.
				pin <- transcripts[which(transcripts[,1]==muh),2]
				pon <- which(row.names(kal)%in%pin)
				pan <- apply(kal[pon,,drop=F], 1, mean)
				# amplicons with TPM > 1
				# get boosted
				# bravo is POSITIVE
				pun <- sapply(pan, function (x)
				{
					if (x>quantile(pan, (as.numeric(PEXPRESSION)/100)))
						{
							as.double((x/max(pan))+100)
						} else { 0 }
				})
				# Here we calc what percentage of the total expression
				# should be captured by each designed pair -> romeo
				# and what score gets calculated according to the
				# significance threshold selected -> bravo
				tango <- do.call(rbind,lapply(strsplit(charlie[,2], '[ _]'), function (f)
				{
					wat <- pun[match(f, names(pan))]
					wit <- pan[match(f, names(pan))]
					wat[which(is.na(wat))] <- 0
					wit[which(is.na(wit))] <- 0
					c(as.numeric(sum(wat)), as.numeric(sum(sapply(wit, function (p)
						{
							(100*p)/sum(pan)
						}))))
				}))
				bravo <- tango[,1]
				romeo <- tango[,2]
				spl['amplicons'] <- psc$detail
				# The eligibility column should report 
				# how many expressed isos are amplified
				spl['expressed_amplified'] <- sapply(1:length(par), function (hi)
				{
					#if (par[hi]!='not_eligible')
					#{
					#	paste0(bravo[hi]%/%100, '.', par[hi], '.', tot)
					#} else {par[hi]}
					if (par[hi]!='not_eligible')
					{
						paste0(bravo[hi]%/%100, '/', sum(pun)%/%100)
					} else {par[hi]}
				})
				november <- as.numeric(lapply(par, function (s)
				{
					if (s=='not_eligible') {-7400} else {eval(parse(text=s))}
				}))
				score <- lapply(1:length(par), function (q)
				{
					bravo[q]-(as.numeric(charlie[q,1])/10)-as.numeric(charlie[q,3])+november[q]
				})
				score <- as.numeric(score)
				# Debugging
				spl['exp_percentage'] <- romeo
				spl['score'] <- score
				#spl['pscoverage'] <- psc$coverage
				#spl['par'] <- par
				#spl['too'] <- too
			}
			spl
		}
		swaggerino
	}))
	stopImplicitCluster()
	stopCluster(cl)
	swagger
}))

swag[,which(sapply(swag, class)=='numeric')]<- apply(swag[,which(sapply(swag, class)=='numeric')],
					    2, round, digits=2)

save(swag, file='finished_product.RData')
system('rm *temp')

# Which pairs do you need
# to amplify all above average isos
fin <- do.call(rbind, lapply(unique(swag$name), function (ua)
	{
		me <- swag[swag$name==ua, ]
		me <- me[order(me$score, decreasing=T), ]
		pin <- transcripts[which(transcripts[,1]==ua),2]
		pon <- which(row.names(kal)%in%pin)
		pan <- apply(kal[pon,,drop=F], 1, mean)
		if (length(pan)!=1)
		{
			pun <- sapply(pan, function (x)
			{
				if (x>quantile(pan, (as.numeric(PEXPRESSION)/100)))
			{ as.double((x/max(pan))+100) } else { 0 }})
			mz <- lapply(names(pun[pun!=0]), function(e) {grep(e, me[me$score>0, 5])})
			me[me$score>0, ][unlist(unique(lapply(mz, function (e) {e[1]}))),]
		} else
		{
			me[me$score>0, ][1, ]
		}
	}))
if (any(is.na(fin$name))) {fin<- fin[-which(is.na(fin$name)),]}

if (nrow(fin)!=0)
{
	# Thorough symbol check with GENCODE annotation
	# you need a gencode table for this
	barcoder <- function (string, feature)
	{
		feature_list <- strsplit(string, '; ')[[1]]
		barcode_pos <- grep(feature, feature_list)
		snippet <- sapply(strsplit(feature_list[barcode_pos], ' '), `[`, 2)
		if (length(barcode_pos)==0) {snippet <- 'NA'}
		return(snippet)
	}
	genid <- lapply(gencode$features[which(gencode$type == 'gene')],
		barcoder, feature='gene_id')
	genid <- as.character(vercutter(genid))
	gename <- lapply(gencode$features[which(gencode$type == 'gene')],
		barcoder, feature='gene_name')
	gename <- as.character(gename)
	genhugo <- lapply(gencode$features[which(gencode$type == 'gene')],
		barcoder, feature='hgnc_id')
	genhugo <- as.character(genhugo)
	gentab <- as.data.frame(cbind(genid, gename, genhugo))
	names(gentab) <- c('name', 'symbol', 'HGNC_support')
	gentab$name<- vercutter(gentab$name)
	fin <- left_join(fin, gentab, by='name')
	fin <- fin[, c(1, 8, 9, 2:4, 6, 7, 5)]
	fin$t_name <- sapply(strsplit(fin$t_name, '_'), `[`, 1)
	fin[, 2:3] <- apply(fin[2:3], 2, as.character)

	ayy <- cbind(sapply(fin$symbol, paste0, '_FOR'), fin$forward, sapply(fin$symbol, paste0, '_REV'), fin$reverse)
	fin2 <- do.call(rbind, apply(ayy, 1, function (q) {as.data.frame(rbind(q[1:2], q[3:4]))}))
	names(fin2)<-c('primer_symbol','primer_sequence')
}

#########################
#                       #
#       CHECKPOINT      #
#                       #
#########################

hs1<- createStyle(fgFill = "#C9D7E9", halign = "CENTER", textDecoration = "Bold",
	border = "TopBottomLeftRight", fontColour = "grey20", borderStyle="thin")
bod<- createStyle(halign='LEFT', border = "TopBottomLeftRight", borderStyle="thin")
xl <- createWorkbook()
if (nrow(fin)!=0)
{
	#write.table(fin, 'validation_candidates.txt', row.names=F, sep='\t', quote=F, na="")
	# save the candidate-for-validation primers 
	addWorksheet(xl, 'validation_candidates')
	writeData(xl, sheet='validation_candidates', x=fin,
		headerStyle=hs1,
		borders="all",
		borderStyle="thin")
	addStyle(xl, sheet='validation_candidates', bod, cols=1:ncol(fin), rows=2:(nrow(fin)+1), gridExpand=T)
	width_vec <- apply(fin, 2, function(x) {max(nchar(as.character(x)) + 4, na.rm = TRUE)})
	width_vec_header <- nchar(colnames(fin)) + 4
	max_vec_header <- pmax(width_vec, width_vec_header)
	setColWidths(xl, sheet='validation_candidates', cols=1:ncol(fin), widths=max_vec_header)
	# save the candidates in a ready for order format
	addWorksheet(xl, 'primers_order')
	writeData(xl, sheet='primers_order', x=fin2,
		headerStyle=hs1,
		borders="all",
		borderStyle="thin")
	addStyle(xl, sheet='primers_order', bod, cols=1:ncol(fin2), rows=2:(nrow(fin2)+1), gridExpand=T)
	width_vec <- apply(fin2, 2, function(x) {max(nchar(as.character(x)) + 4, na.rm = TRUE)})
	width_vec_header <- nchar(colnames(fin2)) + 4
	max_vec_header <- pmax(width_vec, width_vec_header)
	setColWidths(xl, sheet='primers_order', cols=1:ncol(fin2), widths=max_vec_header)
	#write.xlsx(fin, 'validation_candidates.xlsx', asTable = T, overwrite = T, showNA=F)
	#write.xlsx(fin2, 'primers_order.xlsx', asTable = T, overwrite = T, showNA=F)
	# Prep a table for a primersearch check of primers vs the genome
	write.table(fin[, c(3,5,6)], 'val_cand_4genomcheck.txt', row.names=F, sep='\t', quote=F, na="")
	# Launch the in-silico PCR vs the genome
	system('nohup bash psearcher_fingen.sh &')
	system('echo [IsoPrimer] $(date) Generating reports and checking cDNA specificity >> IP_Log.out')
}
# save complete version
addWorksheet(xl, 'primer_omnibus')
writeData(xl, sheet='primer_omnibus', x=swag,
	headerStyle=hs1,
	borders="all",
	borderStyle="thin")
addStyle(xl, sheet='primer_omnibus', bod, cols=1:ncol(swag), rows=2:(nrow(swag)+1), gridExpand=T)
width_vec <- apply(swag, 2, function(x) {max(nchar(as.character(x)) + 4, na.rm = TRUE)})
width_vec_header <- nchar(colnames(swag)) + 4
max_vec_header <- pmax(width_vec, width_vec_header)
setColWidths(xl, sheet='primer_omnibus', cols=1:ncol(swag), widths=max_vec_header)
#write.xlsx(swag, 'primer_omnibus.xlsx', asTable = T, overwrite = T)
saveWorkbook(xl, "primers.xlsx", overwrite=T)

#########################
#                       #
#     REPORT MAKING     #
#                       #
#########################

#######################
#                     #
# MULTIALIGN ISOFORMS #
#                     #
#######################

revcomp<- function (str) {as.character(reverseComplement(DNAString(str)))}

# ClustalOmega of isoforms is substituted to originals
shifter <- function (tapener, fw, rv) 
{
	fuz <- as.integer(aregexec(fw, tapener, max.distance = 5))
	rv <- revcomp(rv)
	zuf <- aregexec(rv, tapener, max.distance = 5)
	zuf <- as.integer(zuf)+nchar(rv)-1
	uber <- mapply(`[`, sapply(tapener, strsplit, ''), mapply(`:`, fuz, zuf, SIMPLIFY=F), SIMPLIFY=F)
	uber_tap <- as.character(lapply(uber, paste0, collapse=''))
	if (length(unique(uber_tap))>1) {
		cw <- msa(uber_tap, type='dna', method='ClustalOmega', order='input')
		cw <- as.character(cw)
		hed <- mapply(`[`, sapply(tapener, strsplit, ''), mapply(`:`, 1, fuz-1, SIMPLIFY=F), SIMPLIFY=F)
		teil <- mapply(`[`, sapply(tapener, strsplit, ''), mapply(`:`, zuf+1, sapply(tapener, nchar), SIMPLIFY=F), SIMPLIFY=F)
		for (e in 1:length(tapener)) {tapener[e] <- paste0(c(hed[[e]], cw[e], teil[[e]]), collapse='')}
	}
	lun<-sapply(tapener, nchar, USE.NAMES=F)
	tail<-lun-as.integer(fuz)
	shift<-lun-tail
	shift<-max(shift)-shift
	ut<-lapply(1:length(tapener), function (g) { paste0(paste0(rep(' ', shift[g]), collapse=''), tapener[g]) })
	sl<-sapply(ut, nchar)
	shift<-max(sl)-sl
	lapply(1:length(tapener), function (g) { paste0(ut[g], paste0(rep(' ', shift[g]), collapse='')) })
}

tagger <- function (t)
{
	sz <- nchar(t)
	tp <- paste0(rep('#', sz+6), collapse='')
	sd <- paste0(c('#  ', rep(' ', sz), '  #'), collapse='')
	nm <- paste0('#  ', t, '  #', collapse='')
	paste0(c(tp, sd, nm, sd, tp), collapse='\n')
}

if (nrow(fin)!=0)
{
	sapply(1:nrow(fin), function (h)
	{
		g <- fin[h, ]
		g <- sapply(g, as.character)
		# Gather data from fin's line
		fw <- g[5]
		rv <- g[6]
		ensid <- g[1]
		gene <- g[3]
		# k is the transcript
		k <- g[4]
		amp <- g[9]
		# Process the amplicon detail
		imp <- strsplit(amp, '[ _]')[[1]][1:(str_count(amp, 'ENST')*2) %% 2==1]
		# Retrieve isoform sequences
		tapener <- cln[which(cln$t_name%in%imp),5]
		tapener <- as.character(gsub(':', '', tapener))
		# Align isoforms
		if (length(tapener)==1) {tap_2 <- tapener} else {
			tap_2 <- shifter(tapener, fw, rv)
		}
		names(tap_2) <- cln[which(cln$t_name%in%imp),3]
		names(tapener) <- names(tap_2)
		# Build tracks
		# what iso the couple was drawn on
		woo <- which(names(tap_2)==k)
		# single base vector
		sq <- strsplit(tap_2[[woo]], '')[[1]]
		# what sj is covered by the couple?
		sj <- lapply(tplus[which(tplus$t_name==k),4], function (x) {as.numeric(strsplit(x, ' ')[[1]])})[[1]]
		sp <- tap_2[[woo]]
		# sj pos on the normalized isoform length
		sh <- unlist(aregexec(fw, sp, max.distance=5))
		# for and rev tracks position
		# and blank tracks
		p1 <- seq(sh, length.out=(nchar(fw)))
		p2 <- seq(unlist(aregexec(revcomp(rv), sp, max.distance=5)), length.out=(nchar(rv)))
		sj <- sj+(sh-regexpr(fw, tapener[woo]))
		bk <- rep(' ', length(sq))
		l1 <- bk
		l2 <- bk
		l4 <- bk
		# actual tracks
		l1[p1] <- strsplit(fw, '')[[1]]
		l1[p2] <- strsplit(revcomp(rv), '')[[1]]
		l2[p1] <- rep('>', nchar(fw))
		l2[p2] <- rep('<', nchar(rv))
		sjn <- which(sj%in%p1|sj%in%p2)
		l4[sj[sjn]] <- '^'
		# multialignment sequence length
		chnk <- 50
		# file name
		fyle <- paste0(c('outputs/', gene, '_aln.txt'), collapse='')
		# gene name
		muh <- ensid
		# kallisto
		pin <- transcripts[which(transcripts[,1]==muh),2]
		pon <- which(row.names(kal)%in%pin)
		pan <- apply(kal[pon,,drop=F], 1, mean)
		# isoforms with expression ≥ mean
		pun <- sapply(pan, function (x) {if (x>quantile(pan, (as.numeric(PEXPRESSION)/100))) { as.double((x/max(pan))+100) } else { 0 }})
		# writing
		write(paste(tagger(gene), '\n\n\n',
			'primers FOR/REV', fw, rv, '\n',
			'designed on isoform:', k, '\n',
			'overlapped sj(s):', paste(sjn, collapse=' '), '\n',
			'most expressed isoform(s):', paste(names(pun)[pun!=0], collapse=' '), '\n',
			'PrimerSearch mismatch:', paste0(PMISMATCH, '%'), '\n',
			'Isoform expression threshold:', paste0(PEXPRESSION, '%'), '\n'), file=fyle, append=T)
		write(paste0('\n', tagger('KALLISTO ESTIMATES'), '\n'), file=fyle, append=T)
		write.table(as.data.frame(pan), file=fyle, col.names=F, row.names=T, sep='\t', quote=F, append=T)
		write(paste0('\n', tagger('MATCHED ISOFORMS MULTIALIGNMENT'), '\n'), file=fyle, append=T)
		for (e in 1:(nchar(tap_2[1])%/%50))
		{
			sb <- ((e-1)*(chnk)+1):(e*chnk)
			write(paste0(c(paste0(l1[sb], collapse=''),
					paste0(l2[sb], collapse=''),
					paste0(lapply(1:length(tap_2), function (i)
					{
						paste(paste0(strsplit(tap_2[[i]], '')[[1]][sb], collapse=''), names(tap_2)[i], sep='\t')
					}), collapse='\n'),
					paste0(l4[sb], collapse=''),
					paste0(l4[sb], collapse=''),
					paste0(l4[sb], collapse='')), collapse='\n'),
						file=fyle, append=T)
		}
		write(paste0('\n', tagger('ISOFORMS FASTA'), '\n'), file=fyle, append=T)
		write.table(tapener, file=fyle, col.names=F, row.names=T, sep='\t', quote=F, append=T)
		write(paste0('\n', tagger('PRIMER3 DATA'), '\n'), file=fyle, append=T)
		gn <- paste0(c(gene, fw, rv), collapse='\t')
		Sys.setenv(GN = gn)
		Sys.setenv(TRANSCRIPT = k)
		system('cd outputs && bash fatail.sh')
		write(paste0('\n', paste0(rep('~', 50), collapse=''), '\n'), file=fyle, append=T)
	})
	system('while [ ! -e genomic_mismatch.txt ]; do sleep 5m; done')
}

# Wait until the specificity check of the primers on the genome is done
# before calling the design run done and removing temp files
if (DEBUG) {system('rm nohup.out sj_seqs.tsv genomic_amplimers.txt val_cand_4genomcheck.txt ./*RData ./*primersearch ./[0-9]*_*txt outputs/[0-9]*_*txt')}

system('echo [IsoPrimer] $(date) Primer design complete >> IP_Log.out')
