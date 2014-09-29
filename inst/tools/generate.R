
debug=TRUE
#debug=FALSE

#cleaned version 2
DUD_URL = "http://dud.docking.org/jahn/DUD_LIB_VS_1.0.tar.gz"
DUDE_URL = "http://dude.docking.org/db/subsets/all/all.tar.gz"
DrugBank_URL = "http://www.drugbank.ca/system/downloads/current/structures/all.sdf.zip"

ensureDataFrame <- function(x) {
	if(is.data.frame(x))
		x
	else if(is.vector(x))
		as.data.frame(as.list(x))
	else{
		warning("ensureDataFrame: found an object of type ",class(x)," where data.frame required")
		x
	}

}
standardFeatures <- function(sdfInput) 
	data.frame(propOB(sdfInput),
				  Ncharges=sapply(bonds(sdfInput, type="charge"), length),
				  #atomcountMA(sdfInput, addH=FALSE),  #variable feature set causes problems
				  ensureDataFrame(groups(sdfInput, type="countMA")), 
				  ensureDataFrame(rings(sdfInput, upper=6, type="count", arom=TRUE)))


buildDud <- function(dbName,downloadDir){

	## Download
	system(paste("wget -c ",DUD_URL,"  -P ",downloadDir))
	system(paste("tar xfz ",file.path(downloadDir,"DUD_LIB_VS_1.0.tar.gz")," -C ",downloadDir))

	## Import
	currentDir= getwd()
	setwd(downloadDir)
	dudDF <- data.frame(actives=paste("./actives/", list.files("./actives/", "sdf$"), sep=""), decoys=paste("./decoys/", list.files("./decoys/", "sdf$"), sep=""))
	row.names(dudDF) <- gsub("^.*actives/|_clustered_.*", "", dudDF[,1])
	setwd(currentDir)

	conn = initDb(dbName)

	loadDud = function(i,type,label)
		loadSdf(conn,file.path(downloadDir,dudDF[i,label]),fct=function(sdfset){
				  data.frame(type=rep(type,length(sdfset)),target_name=rownames(dudDF)[i],standardFeatures(sdfset))
			})

	compoundIds=c()
	for(i in seq(along=dudDF[,1])) {
		compoundIds=c(compoundIds,loadDud(i,"active","actives"))
		compoundIds=c(compoundIds,loadDud(i,"decoy","decoys"))
	}

	compoundIds
}
buildDude <- function(dbName,downloadDir){

	## Download
	system(paste("wget -c ",DUDE_URL,"  -P ",downloadDir))
	system(paste("tar xfz ",file.path(downloadDir,"all.tar.gz")," -C ",downloadDir))
	# read.SDFset(read.SDFstr(gzfile("all/ace/actives <- final.sdf.gz")))

	conn = initDb(dbName)
	currentDir= getwd()
	setwd(downloadDir)
	#actives= list.files(pattern="*actives_final.sdf.gz",recursive=TRUE)
	#decoys= list.files(pattern="*decoys_final.sdf.gz",recursive=TRUE)
	files = list.files(pattern="*_final.sdf.gz",recursive=TRUE)
	if(debug) message("found ",length(files)," .sdf.gz files")
	if(debug) print(head(files))

	compoundIds = c()
	for(file in files){

		if(grepl("actives",file)){
			type="active"
			label="actives"
		}else if(grepl("decoys",file)){
			type="decoy"
			label="decoys"
		}else
			stop("DUDE does not seem to be either an active or decoy")

		if(!file.exists(gsub(".gz$","",file))){
			system(paste("gunzip",file)) #unzip file
			file = gsub(".gz$","",file) #remove .gz extension
			message("loading ", file)

			#read gzipped file
			#sdfData= read.SDFset(read.SDFstr(gzfile(file)))
			name = gsub("all[/\\]|[/\\](actives|decoys)_final.sdf.gz$","",file)
			#if(debug) message("read in ",length(sdfData)," ",type," compounds for ",name)
			compoundIds = c(compoundIds,loadSdf(conn,file,fct=function(sdfset){
						#print(sdfid(sdfset))
						#print(cid(sdfset))
															
					  features = standardFeatures(sdfset)
					  #message("got features: ",paste(colnames(features),collapse=","))
					  #print(head(features))
					  data.frame(type=rep(type,length(sdfset)),target_name=name,features  ) }))
		}

	}
	setwd(currentDir)
	if(debug) message("loaded ",length(compoundIds)," compounds")
	compoundIds

}
buildDrugBank <- function(dbName,downloadDir){


	system(paste("wget -c ",DrugBank_URL,"  -P ",downloadDir))
	system(paste("unzip -d ",downloadDir,file.path(downloadDir,"all.sdf.zip")))

	conn = initDb(dbName)
	loadSdf(conn,file.path(downloadDir,"all.sdf"),fct=standardFeatures)
}
DUD <- function(){
	getDbConn("dud.db")
}
DrugBank <- function(){
	getDbConn("drugbank.db")
}
getDbConn <- function(dbName) {
	initDb(system.file(file.path("extdata",dbName),package="ChemmineDrugs",mustWork=TRUE))
}
