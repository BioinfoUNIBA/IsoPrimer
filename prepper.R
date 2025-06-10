# Use the same annotation version
# as the FASTA you're using
# or some trascripts may go missing
# e.g. gencode.v33.annotation.gtf + gencode.v33.lncRNA_transcripts.fa
# ANNOTATION=gtf table
# THREADS=cores you can run the processes with
# TRANSCRIPTOME=the fasta file

gencode <- read.table(Sys.getenv('ANNOTATION'),
	as.is=T, skip=5, header=F, sep='\t')
gencode <- gencode[,-c(6, 8)]
colnames(gencode) <- c('chr', 'source', 'type', 'start', 'stop', 'strand', 'features')

# Barcoding function: extracts features from the annotation's features column
barcoder <- function (string, feature)
{
	feature_list <- strsplit(string, '; ')[[1]]
	barcode_pos <- grep(feature, feature_list)
	snippet <- sapply(strsplit(feature_list[barcode_pos], ' '), `[`, 2)
	if (length(barcode_pos)==0) {snippet <- 'NA'}
	return(snippet)
}
THREADS<- Sys.getenv('THREADS')
exon_genid <- mclapply(gencode$features[which(gencode$type == 'exon')],
	barcoder, feature='gene_id', mc.cores=THREADS)
exon_genid <- unlist(exon_genid, use.names=F)
exon_trid <- mclapply(gencode$features[which(gencode$type == 'exon')],
	barcoder, feature='transcript_id', mc.cores=THREADS)
exon_trid <- unlist(exon_trid, use.names=F)
exons <- cbind(exon_genid, exon_trid, gencode[which(gencode$type == 'exon'), c(1, 4:5)])
colnames(exons) <- c('name', 't_name', 'chr', 'start', 'stop')
exons['length'] <- exons$stop-exons$start+1

tr_genid <- mclapply(gencode$features[which(gencode$type == 'transcript')],
	barcoder, feature='gene_id', mc.cores=THREADS)
tr_genid <- unlist(tr_genid, use.names=F)
tr_trid <- mclapply(gencode$features[which(gencode$type == 'transcript')],
	barcoder, feature='transcript_id', mc.cores=THREADS)
tr_trid <- unlist(tr_trid, use.names=F)
tr_tsl <- mclapply(gencode$features[which(gencode$type == 'transcript')],
	barcoder, feature='transcript_support_level', mc.cores=THREADS)
tr_tsl <- unlist(tr_tsl, use.names=F)
tr_str <- gencode$strand[which(gencode$type == 'transcript')]

dummy <- table(exon_trid)
tr_exc <- dummy[match(tr_trid, names(dummy))]
tr_exc <- as.vector(tr_exc)
transcripts <- cbind(tr_genid, tr_trid, tr_tsl, tr_exc, gencode[which(gencode$type == 'transcript'), c(4:5)], tr_str)
colnames(transcripts) <- c('name', 't_name', 't_support_lvl', 't_exon_count', 'start', 'stop', 'strand')
tr_len <- aggregate(length~t_name, data=exons, function (x) sum(c(x)))
tr_len <- tr_len[match(tr_trid, tr_len$t_name),]
transcripts['length'] <- tr_len$length

# remove transcripts and exons versions
vercutter <- function (ids)
{
	return(sapply(strsplit(as.character(ids), '\\.'), `[`, 1))
}

transcripts[,1:2] <- apply(transcripts[,1:2], 2, vercutter)
exons[,1:2] <- apply(exons[,1:2], 2, vercutter)

# you just need a table with transcript name, strand and sequence
# for the following to work and mark sjs as ':' in the seq

# Add to the targets.txt file a list of the ensembl IDs of the genes 
# that you want to design primers for
silt <- scan('targets.txt', what='character')
silt <- vercutter(silt)
system(paste0('echo [IsoPrimer.prepper.R] $(date) List of target genes: "', as.character(length(silt)), '" loaded >> IP_Log.out'))

TRANSCRIPTOME<- Sys.getenv('TRANSCRIPTOME')
seqs <- readDNAStringSet(TRANSCRIPTOME)
if (length(grep('transcript_id', names(seqs)))==0)
{
	atnames <- sapply(strsplit(names(seqs), '[^[:alnum:]_]+'), `[`, 1)
} else {
	atnames <- gsub('.*transcript_id=([^]\\.]*).*', '\\1', names(seqs))
}
agnames <- transcripts$name[match(atnames, transcripts$t_name)]
names(seqs) <- paste(atnames, agnames, sep=';')
if (!any(agnames%in%transcripts$name)) {stop(paste0('Please provide the transcriptome FASTA file headers as ">', atnames[1], '"'), call.=FALSE)}
seqx<- seqs[!is.na(agnames)]
writeXStringSet(seqx, 'mod_transcriptome.fa')

seqs <- data.frame(agnames, atnames, as.character(seqs, use.names=F))
seqs<- seqs[!is.na(agnames),]
names(seqs) <- c('name', 't_name', 'sjs')
seqs <- right_join(transcripts[,c(2,7:8)], seqs, by='t_name')
seqs <- seqs[, match(c("name", "strand", "t_name", "length", "sjs"), names(seqs))]
seqs[,c(1,2,3,5)] <- apply(seqs[,c(1,2,3,5)], 2, as.character)

cln <- seqs[which(seqs$name %in% silt), ]
if (nrow(cln)==0) {stop(paste('Please provide target genes in the following format:', agnames[1]), call.=FALSE)}

sj_find <- function (tname, strand, seq)
{
	snp <- exons[which(exons$t_name==tname),]
	if (strand=='+')
	{
		len <- snp$length[order(snp$start)]
	}
	else
	{
		len <- snp$length[order(snp$start, decreasing=T)]
	}
	sj <- sapply(1:length(len), function (e) {sum(len[1:e])+1})
	sj <- sj[-length(sj)]
	sst <- strsplit(seq, "")[[1]]
	sst[sj] <- paste0(':', sst[sj])
	return(paste(sst, collapse=''))
}

sjs <- invisible(sapply(1:nrow(cln),
	function (x)
	{
		sj_find(cln[x, 3], cln[x, 2], cln[x, 5])
	}))

cln <- cbind(cln[,1:4], sjs)
cln$sjs <- as.character(cln$sjs)

# this is your CLN
# "name"   "strand" "t_name" "length" "sjs"
# of the isoforms of all difexp lnc
write.table(cln, 'sj_seqs.tsv', row.names=F, sep='\t', quote=F)

# the tplus for IP.R is 
# "name" "t_name" "transcript_sequence"(without colons) "exons_junction"
# of the isoforms
tplus<-cln
tplus['exons_junction']<-sapply(tplus$sjs, function (x)
	{
		colon <- unlist(gregexpr(':', x))
		paste(colon-1:length(colon), collapse=' ')
	})
tplus$sjs<-sapply(tplus$sjs, function (x) { gsub(':', '', x) })
tplus<-tplus[, c(-2, -4)]
names(tplus)[3]<-'transcript_sequence'

save(cln, tplus, gencode, transcripts, exons, file='wsp.RData')

# proceed with IP.R
