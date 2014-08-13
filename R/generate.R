
#cleaned version 2
DUD_URL = "http://dud.docking.org/jahn/DUD_LIB_VS_1.0.tar.gz"
DrugBank_URL = "http://www.drugbank.ca/system/downloads/current/structures/all.sdf.zip"

standardFeatures <- function(sdfInput) 
	data.frame(propOB(sdfInput),
				  Ncharges=sapply(bonds(sdfInput, type="charge"), length),
				  #atomcountMA(sdfInput, addH=FALSE),  #variable feature set causes problems
				  groups(sdfInput, type="countMA"), 
				  rings(sdfInput, upper=6, type="count", arom=TRUE))


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
	initDb(system.file(file.path("extdata",dbName),package="drugs",mustWork=TRUE))
}
