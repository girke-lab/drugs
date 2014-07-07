
library(ChemmineRdata)

test.buildDud <- function(){
	#DEACTIVATED("off")

	message("building DUD database")
	ids=buildDud("data/dud.db","dataSrc/dud")
	message("loaded ",length(ids)," compounds")

	checkTrue(file.exists("dataSrc/dud/DUD_LIB_VS_1.0.tar.gz"))
	checkTrue(file.exists("dataSrc/dud/actives/vegfr2_clustered_3D_MM.sdf"))
	checkTrue(file.exists("dataSrc/dud/decoys/DUD_vegfr2_decoys_ID_pass_MWPass_I_MM.sdf"))
	checkTrue(file.exists("data/dud.db"))

	
	checkEquals(length(ids),106913)

}
test.buildDrugBank <- function(){
	#DEACTIVATED("off")

	message("building Drug Bank database")
	ids = buildDrugBank("data/drugbank.db","dataSrc/drugbank")
	message("loaded ",length(ids)," compounds")

	checkTrue(file.exists("dataSrc/drugbank/all.sdf"))
	checkTrue(file.exists("data/drugbank.db"))

	checkEquals(length(ids),6795) # 6806 - 11 failures

}
