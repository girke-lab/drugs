# 7109 7087
library(ChemmineDrugs)


targetFeatures = sort(c(
	"aromatic","cansmi","cansmins","formula",
	"hba1","hba2", "hbd","inchi","logp","mr",
	"mw","ncharges", "nf","r2nh","r3n","rcch",
	"rcho","rcn", "rcooh","rcoor","rcor","rings",
	"rnh2","roh", "ropo3","ror","title","tpsa"
))

test.buildDude <- function(){
	DEACTIVATED("does not work as a test anymore. only here for reference")

	message("building DUDE database")
	ids=ChemmineDrugs:::buildDude("inst/extdata/dud.db","dataSrc/dud")
	message("loaded ",length(ids)," compounds")

	checkTrue(file.exists("dataSrc/dud/all.tar.gz"))
	checkTrue(file.exists("dataSrc/dud/all/ace/actives_final.sdf.gz"))
	checkTrue(file.exists("dataSrc/dud/all/plk1/decoys_final.sdf.gz"))
	checkTrue(file.exists("inst/extdata/dud.db"))

	
	# this will only succeed on a clean load
	#checkEquals(length(ids),1449671)
	#1454636 file
	#1449671 db

	conn = initDb("inst/extdata/dud.db")
	allIds = getAllCompoundIds(conn)
	message("found ",length(allIds)," compounds")
	checkEquals(length(allIds),1449670)


}

test.buildDrugBank <- function(){
	DEACTIVATED("does not work as a test anymore. only here for reference")

	message("building Drug Bank database")
	ids = ChemmineDrugs:::buildDrugBank("inst/extdata/drugbank.db","dataSrc/drugbank")
	message("loaded ",length(ids)," compounds")

	checkTrue(file.exists("dataSrc/drugbank/all.sdf"))
	checkTrue(file.exists("inst/extdata/drugbank.db"))

	checkEquals(length(ids),6795) # 6806 - 11 failures

}
test.zzDudFeatures<- function(){
	#DEACTIVATED("off")
	features = listFeatures(DUD())
	print(features)
	checkEquals(sort(features),sort(c(targetFeatures,"type","target_name")))
}
test.zzDrugBankFeatures<- function(){
	#DEACTIVATED("off")
	features = listFeatures(DrugBank())
	print(features)
	checkEquals(sort(features),targetFeatures)
}
