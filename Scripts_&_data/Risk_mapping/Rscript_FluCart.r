library(blockCV)
library(dismo)
library(gbm)
library(lubridate)
library(raster)
library(RColorBrewer)
library(sf)

# 1. Generating a raster reflecting the density of close wild bird cases
# 2. Loading and preparing the environmental layers and occurrence data
# 3. Ecological niche modelling analyses based on presence/absence data
# 4. Ecological niche modelling analyses based on the number of herds
# 5. Analysing the ecological niche models (RI, response curves, mapping)
# 6. Generating a composite figure with the different risk maps

savingFiles = FALSE; savingPlots = FALSE; nberOfReplicates = 100
dataset = read.csv("Data_160625.csv", head=T, sep=";"); dates = dataset[,"final_date"]
dates = decimal_date(dmy(dates)); minYear = min(dates); maxYear = max(dates)
analyses = c("a","b","c") # based on farm and hobby (a), farm (b), and hobby (c) data
analyse_names = c("Contamined farm and hobby units", "Contamined farms", "Contamined hobby units")

# 1. Generating a raster reflecting the density of close wild bird cases
	 # (positive wild bird cases per km2 in the municipality and adjacent ones)

emptyRaster = raster("All_env_variables/Accessibility/TTC12_res_BELG_ETRS_LAEA.tif")
emptyRaster[which(!is.na(emptyRaster[]))] = 0; contour = rasterToPolygons(emptyRaster>0, dissolve=T)
if (!file.exists("All_env_variables/Wild_bird_cases_density/Wild_bird_cases_density.tif"))
	{
		municipalities = shapefile("BE_municipalities/Shapefile_NIS5_codes.shp")
		municipalities_wgs = spTransform(municipalities, CRS("+init=epsg:4326"))
		municipalities_lam = spTransform(municipalities, crs(emptyRaster))
		tab = data.frame(matrix(nrow=dim(municipalities@data)[1], ncol=5))
		colnames(tab) = c("NIS5","area","adjacent_municipalities","number_of_wild_bird_cases",
						  "wild_bird_cases_density_in_all_adjacent_municipalities")
		tab[,"NIS5"] = municipalities_lam@data[,"NISCode"]
		tab[,"area"] = raster::area(municipalities_wgs)/1000000 # in km2
			# sum(tab[,"area"]) = 30689.43 km2 (which is correct for Belgium)
		intersections = st_intersects(st_as_sf(municipalities_lam), st_as_sf(municipalities_lam))
		for (i in 1:dim(municipalities_lam@data)[1])
			{
				indices = which((dataset[,"animal_type"]=="Wild")
								&(dataset[,"NIS_code"]==municipalities_lam@data[i,"NISCode"])
								&(dataset[,"Result_AFLoccurence"]=="POSITIVE"))
				tab[i,"number_of_wild_bird_cases"] = length(indices); buffer = intersections[[i]][1]
				for (j in 2:length(intersections[[i]])) buffer = paste0(buffer,"-",intersections[[i]][j])
				tab[i,"adjacent_municipalities"] = buffer
			}
		for (i in 1:dim(municipalities@data)[1])
			{
				total_cases = as.numeric(tab[i,"number_of_wild_bird_cases"])
				total_area = as.numeric(tab[i,"area"])
				indices = as.numeric(unlist(strsplit(tab[i,"adjacent_municipalities"],"-")))
				for (j in 1:length(indices))
					{
						total_cases = total_cases + tab[indices[j],"number_of_wild_bird_cases"]
						total_area = total_area + tab[indices[j],"area"]
					}
				tab[i,"wild_bird_cases_density_in_all_adjacent_municipalities"] = total_cases/total_area
			}
		write.csv(tab, "BE_municipalities/Wild_bird_cases_NIS5.csv", row.names=F, quote=F)
		tab = read.csv("BE_municipalities/Wild_bird_cases_NIS5.csv", head=T); buffer = emptyRaster
		points = data.frame(xyFromCell(emptyRaster,1:length(emptyRaster[]))); colnames(points) = c("x","y")
		indices = st_intersects(st_as_sf(points, coords=c("x","y"), crs=crs(municipalities_lam)), st_as_sf(municipalities_lam))
		for (i in 1:length(indices))
			{
				if (length(indices[[i]]) != 0) buffer[i] = tab[indices[[i]],"wild_bird_cases_density_in_all_adjacent_municipalities"]
			}
		# plot(buffer); points(coordinates(municipalities_lam), cex=tab[,"wild_bird_cases_density_in_all_adjacent_municipalities"]*10)
		writeRaster(buffer, "All_env_variables/Wild_bird_cases_density/Wild_bird_cases_density.tif")	
	}
if (!file.exists("All_env_variables/Farms_&_poultry_2020-23/Sciensano_farms_density.tif"))
	{
		rast1 = raster("All_env_variables/Farms_&_poultry_2020-23/Sciensano_farms_density_1x1k.tif")
		rast2 = raster::aggregate(rast1, 5, fun=sum)
		rast3 = raster::projectRaster(rast2, emptyRaster)
		rast3[which(!is.na(emptyRaster[]))] = round(rast3[which(!is.na(emptyRaster[]))])
		rast3[which(is.na(emptyRaster[]))] = NA
		writeRaster(rast3, "All_env_variables/Farms_&_poultry_2020-23/Sciensano_farms_density.tif", overwrite=T)
		rast1 = raster("All_env_variables/Farms_&_poultry_2020-23/Sciensano_poultry_population_1x1k.tif")
		rast2 = raster::aggregate(rast1, 5, fun=sum)
		rast3 = raster::projectRaster(rast2, emptyRaster)
		rast3[which(!is.na(emptyRaster[]))] = round(rast3[which(!is.na(emptyRaster[]))])
		rast3[which(is.na(emptyRaster[]))] = NA
		writeRaster(rast3, "All_env_variables/Farms_&_poultry_2020-23/Sciensano_poultry_population.tif", overwrite=T)
	}

# 2. Loading and preparing the environmental layers and occurrence data

	# 2.1. Loading and preparing the environmental rasters

envVariables = list(); envVariableNames1 = list(); envVariableNames2 = list()
i = 1; envVariables[[i]] = raster("All_env_variables/Accessibility/TTC12_res_Belg_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Inaccessibility"
i = 2; envVariables[[i]] = raster("All_env_variables/Distance_to_coast/DTC_AsCrowFlies_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Distance to coast (km)"
i = 3; envVariables[[i]] = raster("All_env_variables/Human_population_density/Pop_density_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Human population"
i = 4; envVariables[[i]] = raster("All_env_variables/Farms_&_poultry_2020-23/Sciensano_farms_density.tif"); envVariableNames2[[i]] = "Farms density"
i = 5; envVariables[[i]] = raster("All_env_variables/Farms_&_poultry_2020-23/Sciensano_poultry_population.tif"); envVariableNames2[[i]] = "Poultry population"
i = 6; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_Bfor_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Broadleaf forests (%)"
i = 7; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_Cfor_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Coniferous forests (%)"
i = 8; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_Crop_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Croplands (%)"
i = 9; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_for_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "All forest types (%)"
i = 10; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_Grass_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Grasslands (%)"
i = 11; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_Mfor_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Mixed forests (%)"
i = 12; envVariables[[i]] = raster("All_env_variables/Land_cover/Perc_Urb_ref_grid_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Urban areas (%)"
i = 13; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Air_temperature/Fall_Tmean_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Avg. temperature (fall)"
i = 14; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Air_temperature/Spring_Tmean_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Avg. temperature (spring)"
i = 15; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Air_temperature/Summer_Tmean_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Avg. temperature (summer)"
i = 16; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Air_temperature/Winter_Tmean_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Avg. temperature (winter)"
i = 17; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Cumulated_rainfall/Fall_Crain_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Cum. precipitation (fall)"
i = 18; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Cumulated_rainfall/Spring_Crain_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Cum. precipitation (spring)"
i = 19; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Cumulated_rainfall/Summer_Crain_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Cum. precipitation (summer)"
i = 20; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Cumulated_rainfall/Winter_Crain_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Cum. precipitation (winter)"
i = 21; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Nb_rainy_days/Fall_Nb_RD_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Rainfall days (fall)"
i = 22; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Nb_rainy_days/Spring_Nb_RD_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Rainfall days (spring)"
i = 23; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Nb_rainy_days/Summer_Nb_RD_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Rainfall days (summer)"
i = 24; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Nb_rainy_days/Winter_Nb_RD_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Rainfall days (winter)"
i = 25; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Relative_humidity/Fall_RH_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Average RH (fall)"
i = 26; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Relative_humidity/Spring_RH_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Average RH (spring)"
i = 27; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Relative_humidity/Summer_RH_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Average RH (summer)"
i = 28; envVariables[[i]] = raster("All_env_variables/Meteo_data_2021-23/Relative_humidity/Winter_RH_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Average RH (winter)"
i = 29; envVariables[[i]] = raster("All_env_variables/Water_bodies/Top10Vector_Perc_WC_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Water courses"
i = 30; envVariables[[i]] = raster("All_env_variables/Water_bodies/Top10Vector_Perc_WL_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Wetlands"
i = 31; envVariables[[i]] = raster("All_env_variables/Water_bodies/Top10Vector_Perc_WS_ETRS_LAEA.tif"); envVariableNames2[[i]] = "Water surfaces"
i = 32; envVariables[[i]] = raster("All_env_variables/Wild_bird_cases_density/Wild_bird_cases_density.tif"); envVariableNames2[[i]] = "Density of close wild bird cases"
for (i in 1:length(envVariables)) envVariableNames1[[i]] = names(envVariables[[i]])
i = 4; envVariableNames1[[i]] = "Sciensano_farms_density"; names(envVariables[[i]]) = envVariableNames1[[i]]
i = 5; envVariableNames1[[i]] = "Sciensano_poultry_population"; names(envVariables[[i]]) = envVariableNames1[[i]]
nullRaster = envVariables[[1]]; nullRaster[which(!is.na(nullRaster[]))] = 1
envVariableNames1 = unlist(envVariableNames1); envVariableNames2 = unlist(envVariableNames2)

	# 2.2. Loading and preparing the occurrence data

coordinates_lam = matrix(nrow=dim(dataset)[1], ncol=2); colnames(coordinates_lam) = c("longitude_Lambert","latitude_Lambert")
for (i in 1:dim(dataset)[1])
	{
		if (!is.na(dataset[i,"longitude_all"]))
			{
				buffer = data.frame(lon=as.numeric(dataset[i,"longitude_all"]), lat=as.numeric(dataset[i,"latitude_all"]))
				coordinates(buffer) = c("lon", "lat"); crs(buffer) = CRS("+init=epsg:4326") # WGS84
				coordinates_lam[i,c("longitude_Lambert","latitude_Lambert")] = spTransform(buffer, crs(nullRaster))@coords
			}
	}
if (!"longitude_Lambert"%in%colnames(dataset)) dataset = cbind(dataset, coordinates_lam)
temp = dataset[which(is.na(raster::extract(nullRaster,dataset[,c("longitude_Lambert","latitude_Lambert")]))),]
if (savingFiles) write.csv(temp, "Outside_BEL.csv", row.names=F, quote=F)
dataset = dataset[which(!is.na(raster::extract(nullRaster,dataset[,c("longitude_Lambert","latitude_Lambert")]))),]

# 3. Ecological niche modelling analyses based on presence/absence data
	
	# 3.1. Generating the BRT dataframes (one per replicate analysis)

envVariablesWithNotEnoughVariation = c("Perc_Bfor_ref_grid_ETRS_LAEA","Perc_Cfor_ref_grid_ETRS_LAEA")
if (!file.exists(paste0("Dataframes_1.rds")))
	{
		dataframes_list = list()
		for (h in 1:length(analyses))
			{
				if (analyses[h] == "a") data = dataset[which(dataset[,"animal_type"]%in%c("Farm","Hobby","Hobby?")),]
				if (analyses[h] == "b") data = dataset[which(dataset[,"animal_type"]%in%c("Farm")),]
				if (analyses[h] == "c") data = dataset[which(dataset[,"animal_type"]%in%c("Hobby","Hobby?")),]
				positive_entries = which((data$Result_H5H7Other=="H5")&(data$Result_Pathotype=="HP")&(data$Result_AFLoccurence=="POSITIVE"))
				positive_entries_coordinates = data[positive_entries,c("longitude_Lambert","latitude_Lambert")]
				positive_entries_coordinates_text = paste0(positive_entries_coordinates[,1],"_",positive_entries_coordinates[,2])
				positive_entries_coordinates_text_unique = unique(positive_entries_coordinates_text)
				positive_entries_coordinates_unique = matrix(nrow=length(positive_entries_coordinates_text_unique), ncol=2)
				for (i in 1:length(positive_entries_coordinates_text_unique))
					{
						positive_entries_coordinates_unique[i,] = unlist(strsplit(positive_entries_coordinates_text_unique[i],"_"))
					}
				positive_entries_coordinates_unique = positive_entries_coordinates_unique[which(positive_entries_coordinates_unique[,1]!="NA"),]
				positive_entries_coordinates_unique = as.data.frame(positive_entries_coordinates_unique)
				positive_entries_coordinates_unique[,1] = as.numeric(positive_entries_coordinates_unique[,1])
				positive_entries_coordinates_unique[,2] = as.numeric(positive_entries_coordinates_unique[,2])
				positive_cells_unique = unique(cellFromXY(nullRaster, positive_entries_coordinates_unique))
				positive_entries_coordinates_cells = xyFromCell(nullRaster, positive_cells_unique)
				negative_entries = which(data$Result_AFLoccurence=="NEGATIVE")
				negative_entries_coordinates = data[negative_entries,c("longitude_Lambert","latitude_Lambert")]
				negative_entries_coordinates_text = paste0(negative_entries_coordinates[,1],"_",negative_entries_coordinates[,2])
				negative_entries_coordinates_text_unique = unique(negative_entries_coordinates_text)
				negative_entries_coordinates_unique = matrix(nrow=length(negative_entries_coordinates_text_unique), ncol=2)
				for (i in 1:length(negative_entries_coordinates_text_unique))
					{
						negative_entries_coordinates_unique[i,] = unlist(strsplit(negative_entries_coordinates_text_unique[i],"_"))
					}
				negative_entries_coordinates_unique = negative_entries_coordinates_unique[which(negative_entries_coordinates_unique[,1]!="NA"),]
				negative_entries_coordinates_unique = as.data.frame(negative_entries_coordinates_unique)
				negative_entries_coordinates_unique[,1] = as.numeric(negative_entries_coordinates_unique[,1])
				negative_entries_coordinates_unique[,2] = as.numeric(negative_entries_coordinates_unique[,2])
				negative_cells_unique = unique(cellFromXY(nullRaster, negative_entries_coordinates_unique))
				negative_cells_unique = negative_cells_unique[which(!negative_cells_unique%in%positive_cells_unique)]
				negative_entries_coordinates_cells = xyFromCell(nullRaster, negative_cells_unique)
				df1a = cbind(rep(1,dim(positive_entries_coordinates_cells)[1]),positive_entries_coordinates_cells)
				colnames(df1a) = c("response","x","y"); df1b = as.data.frame(df1a); dataframes = list()
				for (i in 1:length(envVariables))
					{
						if (!envVariableNames1[i]%in%envVariablesWithNotEnoughVariation)
							{
								df1b[,names(envVariables[[i]])] = raster::extract(envVariables[[i]], df1a[,c("x","y")])
							}
					}
				for (i in 1:nberOfReplicates)
					{
						df2a = df1a; df2a[,"response"] = 0; df2a[,c("x","y")] = NA; indices = 1:dim(negative_entries_coordinates_cells)[1]
						df2a[,c("x","y")] = negative_entries_coordinates_cells[sample(indices,dim(df1a)[1],replace=F),c("x","y")]
						df2b = as.data.frame(df2a)
						for (j in 1:length(envVariables))
							{
								if (!envVariableNames1[j]%in%envVariablesWithNotEnoughVariation)
									{
										df2b[,names(envVariables[[j]])] = raster::extract(envVariables[[j]], df2a[,c("x","y")])
									}
							}
						dataframes[[i]] = rbind(df1b, df2b)
					}
				for (i in 1:nberOfReplicates)
					{
						for (j in 4:dim(dataframes[[i]])[2])
							{
								if (length(unique(dataframes[[i]][,j])) == 1) print(c(h,i,j-3))
							}
					}
				dataframes_list[[h]] = dataframes
			}
		saveRDS(dataframes_list, paste0("Dataframes_1.rds"))
	}	else	{
		dataframes_list = readRDS(paste0("Dataframes_1.rds"))
	}

	# 3.2. Investigating the optimal spatial fold size (with "blockCV")

if (!file.exists(paste0("Correlograms.pdf")))
	{
		pdf(paste0("Correlograms_NEW.pdf"), width=4.5, height=4); par(mfrow=c(3,1), mar=c(2,2,1,1), oma=c(0.5,0.5,0,0))
		for (h in 1:length(analyses))
			{
				dataframes = dataframes_list[[h]]; correlograms = list(); values = rep(NA, nberOfReplicates)
				for (i in 1:nberOfReplicates)
					{
						correlograms[[i]] = ncf::correlog(dataframes[[i]][,"x"], dataframes[[i]][,"y"], dataframes[[i]][,"response"], na.rm=T, increment=5000, resamp=0, latlon=F)
					}
				plot(correlograms[[1]]$mean.of.class, correlograms[[1]]$correlation, ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, xlim=c(3000,200000), ylim=c(-0.75,0.75))
				for (i in 1:nberOfReplicates)
					{
						values[i] = correlograms[[i]]$mean.of.class[which(correlograms[[i]]$correlation <= 0)[1]]
						lines(correlograms[[i]]$mean.of.class, correlograms[[i]]$correlation, lwd=0.1, col="gray60")
					}
				meanV = round(mean(values)); IC95 = round(quantile(values, c(0.025,0.975)))
				cat(meanV,", 2.5th-97.5th percentiles range = [",IC95[1],",",IC95[2],"]\n",sep=""); abline(h=0, lwd=0.5, col="red", lty=2)
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,250000,50000))
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-1,1,0.5))
				title(xlab="Distance (m)", cex.lab=0.65, mgp=c(0.9,0,0), col.lab="gray30")
				title(ylab="Correlation", cex.lab=0.65, mgp=c(1.2,0,0), col.lab="gray30")
				mtext(toupper(analyses[h]), side=3, line=-0.7, at=30000, cex=0.6, col="gray30")
			}
		dev.off() # no obvious spatial autocorrelation --> arbitrary selection of a distance of ~30000 (30km)
	}
ranges = c(25000, 30000, 10000)
blockSizeExplorationByBlockCV = FALSE
if (blockSizeExplorationByBlockCV == TRUE)
	{
		buffer = list(); ranges1 = matrix(nrow=length(analyses), ncol=1); row.names(ranges1) = analyses
		for (i in 1:length(envVariables)) buffer[[i]] = terra::rast(envVariables[[i]])
		for (i in 1:length(envVariables)) names(buffer[[i]]) = envVariableNames1[i]
		test = cv_spatial_autocor(r=terra::rast(buffer), # a SpatRaster object or path to files
								  num_sample=1000, plot=T) # num_sample: number of cells to be used
		for (i in 1:length(analyses))
			{
				ranges2 = rep(NA, nberOfReplicates) 
				for (j in 1:nberOfReplicates)
					{
						pa_data = sf::st_as_sf(data.frame(dataframes_list[[i]][[j]]), coords=c("x","y"), crs=crs(emptyRaster))
						trycatch = tryCatch(
							{
								ranges2[j] = cv_spatial_autocor(x=pa_data, column="response", plot=F)$range
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
					}
				cat(i,", ",gsub("_","",analyses[i]),": ",length(which(is.na(ranges2)))," \"NA\"","\n",sep="")
				ranges1[i,1] = round(mean(ranges2, na.rm=T))
			}
		ranges = ranges1
	}

	# 3.3. Training all the BRT models with two cross-validations
	
newAnalyses = FALSE
keep.fold.vector = FALSE; keep.fold.fit = FALSE; subdirectory = "BRT_models"
# keep.fold.vector = TRUE; keep.fold.fit = TRUE; subdirectory = "BRT_foldFit"
if (newAnalyses == TRUE) { for (h in 1:length(analyses))	 {
		rasters_stack = stack(envVariables)
		names(rasters_stack) = envVariableNames1
		theRanges = ranges[h] # distance in meters
		gbm.x = gsub("\\.","-",names(rasters_stack))
		gbm.x = gbm.x[which(!gbm.x%in%envVariablesWithNotEnoughVariation)]
		gbm.y = "response"
		offset = NULL
		tree.complexity = 5 # "tc" = number of nodes in the trees
		learning.rate = 0.001 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(dataframes_list[[h]][[1]])[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "bernoulli"
		n.trees = 10 # initial number of trees
		step.size = 2 # interval at which the predictive deviance is computed and logged
					  # (at each interval, the folds are successively used as test data set
					  # nd the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		showingFoldsPlot = FALSE
		brt_models_ccv = list() # classic cross-validations
		brt_models_scv = list() # spatial cross-validations
		AUCs = matrix(nrow=nberOfReplicates, ncol=2); colnames(AUCs) = c("CCV_AUC","SCV_AUC")
		for (i in 1:nberOfReplicates)
			{
				data = as.data.frame(dataframes_list[[h]][[i]]); n.trees = 10; learning.rate = 0.001; step.size = 2; fold.vector = NULL; worked = FALSE
				# BRT with classic (standard) cross-validation (CCV):
				pdf(file=paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_model_",analyses[h],"_CCV_",i,".pdf"))
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								brt_models_ccv[[i]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
						if (length(brt_models_ccv) == i) worked = TRUE
					}
				dev.off()
				AUCs[i,"CCV_AUC"] = brt_models_ccv[[i]]$cv.statistics$discrimination.mean
				# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
				spdf = SpatialPointsDataFrame(data[c("x","y")], data[,c(1,4:dim(data)[2])], proj4string=crs(nullRaster)); worked = FALSE
				pdf(file=paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_model_",analyses[h],"_SCV_",i,".pdf")); worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								myblocks = cv_spatial(spdf, column="response", k=n.folds, size=theRanges[1], selection="random")
								n.trees = 10; learning.rate = 0.001; step.size = 2; fold.vector = myblocks$folds_ids # "myblocks$foldID" not correct!
								brt_models_scv[[i]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
						if (length(brt_models_scv) == i) worked = TRUE
					}
				dev.off()
				AUCs[i,"SCV_AUC"] = brt_models_scv[[i]]$cv.statistics$discrimination.mean
			}
		saveRDS(brt_models_ccv, paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_models_",analyses[h],"_CCV.rds"))
		saveRDS(brt_models_scv, paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds"))		
		write.csv(AUCs, paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_models_",analyses[h],"_AUC.csv"), row.names=F, quote=F)
	}}

	# 3.4. Computation of the Sørensen indexes (SIs) and associated optimsed thresholds

SIs = matrix(nrow=nberOfReplicates, ncol=6); tabs_list = list(); subdirectory = "BRT_foldFit"
colnames(SIs) = c(paste0(rep("SI_",3),analyses), paste0(rep("threshold_",3),analyses))
for (h in 1:length(analyses))
	{
		brt_models_scv = readRDS(paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds")); tabs = list()
		for (i in 1:length(brt_models_scv))
			{
				responses = brt_models_scv[[i]]$gbm.call$dataframe[,gbm.y]; foldFit = brt_models_scv[[i]]$fold.fit
				prediction = stats::plogis(foldFit) # back-transformation from the logit to the probability scale
				P = sum(responses==1); A = sum(responses==0); prev = P/(P+A); x = (P/A)*((1-prev)/prev)
				SI = 0; optimised_threshold = NA; tmp = matrix(nrow=101, ncol=2); tmp[,1] = seq(0,1,0.01)
				for (threshold in seq(0,1,0.01))
					{
						TP = length(which((responses==1)&(prediction>=threshold))) # true positives
						FN = length(which((responses==1)&(prediction<threshold))) # false negatives
						FP_pa = length(which((responses==0)&(prediction>=threshold))) # false positives
						SI_tmp = (2*TP)/((2*TP)+(x*FP_pa)+(FN))
						if (SI < SI_tmp)
							{
								SI = SI_tmp; optimised_threshold = threshold
							}
						tmp[which(tmp[,1]==threshold),2] = SI_tmp
					}
				SIs[i,paste0("SI_",analyses[h])] = round(SI,3)
				SIs[i,paste0("threshold_",analyses[h])] = optimised_threshold
				tabs[[i]] = tmp
			}
		tabs_list[[h]] = tabs
		for (i in 1:dim(SIs)[1])
			{
				if (nchar(SIs[i,paste0("SI_",analyses[h])]) == 3) SIs[i,paste0("SI_",analyses[h])] = paste0(SIs[i,paste0("SI_",analyses[h])],"00")
				if (nchar(SIs[i,paste0("SI_",analyses[h])]) == 4) SIs[i,paste0("SI_",analyses[h])] = paste0(SIs[i,paste0("SI_",analyses[h])],"0")
				if (!grepl("\\.",SIs[i,paste0("threshold_",analyses[h])])) SIs[i,paste0("threshold_",analyses[h])] = paste0(SIs[i,paste0("threshold_",analyses[h])],".00")
				if (nchar(SIs[i,paste0("threshold_",analyses[h])]) != 4) SIs[i,paste0("threshold_",analyses[h])] = paste0(SIs[i,paste0("threshold_",analyses[h])],"0")
			}
		write.csv(SIs, "SI_thresholds.csv", row.names=F, quote=F); vS = as.numeric(SIs[,paste0("SI_",analyses[h])])
		quantiles = round(quantile(vS,c(0.025,0.975)),3); meanV = round(mean(vS),3)
		cat("\tanalysis (",analyses[h],"): mean SI = ",meanV," (2.5th-97.5th percentiles range = [",quantiles[1],"-",quantiles[2],"])\n",sep="")
			# analysis (a): mean SI = 0.676 (2.5th-97.5th percentiles range = [0.667-0.733])
			# analysis (b): mean SI = 0.710 (2.5th-97.5th percentiles range = [0.667-0.779])
			# analysis (c): mean SI = 0.675 (2.5th-97.5th percentiles range = [0.667-0.726])
	}
if (savingPlots)
	{
		pdf(paste0("All_SI_curves.pdf"), width=6.3, height=1.95)
		par(mfrow=c(1,3), oma=c(0,0.3,0,0), mar=c(3,3,0.5,0.5), lwd=0.4, col="gray30")
		for (h in 1:length(analyses))
			{
				plot(tabs_list[[h]][[1]], col=NA, ann=F, axes=F, xlim=c(0,1), ylim=c(0,1))
				for (i in 1:length(tabs_list[[h]])) lines(tabs_list[[h]][[i]], lwd=0.3, col="gray80", lty=1)
				axis(side=1, lwd.tick=0.4, cex.axis=0.8, lwd=0, tck=-0.03, col.axis="gray30", col.ticks="gray30", mgp=c(0,0.25,0))
				axis(side=2, lwd.tick=0.4, cex.axis=0.8, lwd=0, tck=-0.03, col.axis="gray30", col.ticks="gray30", mgp=c(0,0.45,0))
				title(ylab=expression("SI"), cex.lab=1, mgp=c(1.5,0,0), col.lab="gray30")
				title(xlab="Threshold", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30"); box(lwd=0.4, col="gray30")
				mtext(toupper(analyses[h]), side=3, line=-1.40, at=0.04, cex=0.7, col="gray30")
			}
		dev.off()
	}

# 4. Ecological niche modelling analyses based on the number of herds

	# 4.1. Generating the BRT dataframes (one per replicate analysis)

envVariablesWithNotEnoughVariation = c("Perc_Bfor_ref_grid_ETRS_LAEA","Perc_Cfor_ref_grid_ETRS_LAEA","Perc_WB_ref_grid_ETRS_LAEA")
if (!file.exists(paste0("Dataframes_2.rds")))
	{
		dataframes_list = list()
		for (h in 1:length(analyses))
			{
				if (analyses[h] == "a") data = dataset[which(dataset[,"animal_type"]%in%c("Farm","Hobby","Hobby?")),]
				if (analyses[h] == "b") data = dataset[which(dataset[,"animal_type"]%in%c("Farm")),]
				if (analyses[h] == "c") data = dataset[which(dataset[,"animal_type"]%in%c("Hobby","Hobby?")),]
				unit_IDs = data[,"animal_id_group2"]
				for (i in 1:dim(data)[1])
					{
						if (data[i,"animal_type"] == "Farm") unit_IDs[i] = data[i,"farm_id_all_uniq"]
					}
				unit_IDs = unique(unit_IDs); unit_IDs = unit_IDs[order(unit_IDs)]
				units = matrix(nrow=length(unit_IDs), ncol=7); units[,1] = unit_IDs
				colnames(units) = c("unit_ID","cell_ID","x","y","unit_type","test","pathotype")
				units = as.data.frame(units)
				for (i in 1:length(unit_IDs))
					{
						index = which(data[,"animal_id_group2"]==unit_IDs[i])
						if (length(index) > 0)
							{
								index = index[1]
							}	else	{
								index = which(data[,"farm_id_all_uniq"]==unit_IDs[i])[1]
							}
						units[i,"cell_ID"] = cellFromXY(nullRaster, data[index,c("longitude_Lambert","latitude_Lambert")])
						units[i,c("x","y")] = xyFromCell(nullRaster, as.numeric(units[i,"cell_ID"]))
						units[i,"unit_type"] = stringr::str_to_lower(data[index,"animal_type"])
						units[i,"test"] = stringr::str_to_lower(data[index,"Result_AFLoccurence"])
						units[i,"pathotype"] = data[index,"Result_H5H7Other"]
					}
				units = units[which(!grepl("H7",units[,"pathotype"])),]
				cell_IDs = as.numeric(unique(units[,"cell_ID"])); cell_IDs = cell_IDs[order(cell_IDs)]
				cells = matrix(nrow=length(cell_IDs), ncol=4)
				colnames(cells) = c("response","cell_ID","x","y")
				cells = as.data.frame(cells); cells[,"cell_ID"] = cell_IDs
				for (i in 1:length(cell_IDs))
					{
						index = which(units[,"cell_ID"]==cell_IDs[i])[1]
						cells[i,"x"] = as.numeric(units[index,"x"])
						cells[i,"y"] = as.numeric(units[index,"y"])
						cells[i,"response"] = sum((units[,"test"]=="positive")&(units[,"cell_ID"]==cell_IDs[i]))
					}
				df1a = cells[which(cells[,"response"]>0),c("response","x","y")]
				df1a = as.data.frame(df1a); df1b = df1a; dataframes = list()
				for (i in 1:length(envVariables))
					{
						if (!envVariableNames1[i]%in%envVariablesWithNotEnoughVariation)
							{
								df1b[,names(envVariables[[i]])] = raster::extract(envVariables[[i]], df1a[,c("x","y")])
							}
					}
				for (i in 1:nberOfReplicates)
					{
						buffer = as.data.frame(units[which(units[,"test"]=="negative"),])
						buffer$x = as.numeric(buffer$x); buffer$y = as.numeric(buffer$y)
						df2a = df1a; df2a[,"response"] = 0; df2a[,c("x","y")] = NA; indices = 1:dim(buffer)[1]
						df2a[,c("x","y")] = buffer[sample(indices,dim(df1a)[1],replace=F),c("x","y")]; df2b = as.data.frame(df2a)
						for (j in 1:length(envVariables))
							{
								if (!envVariableNames1[j]%in%envVariablesWithNotEnoughVariation)
									{
										df2b[,names(envVariables[[j]])] = raster::extract(envVariables[[j]], df2a[,c("x","y")])
									}
							}
						dataframes[[i]] = rbind(df1b, df2b)
					}
				dataframes_list[[h]] = dataframes
			}
		saveRDS(dataframes_list, "Dataframes_2.rds")
	}	else	{
		dataframes_list = readRDS("Dataframes_2.rds")
	}

	# 4.2. Training all the BRT models with two cross-validations
	
keep.fold.vector = FALSE; keep.fold.fit = FALSE; subdirectory = "BRT_models"
# keep.fold.vector = TRUE; keep.fold.fit = TRUE; subdirectory = "BRT_foldFit"
newAnalyses = FALSE # note: SCV not working for analyses #2 and #3
if (newAnalyses == TRUE) { for (h in 1:length(analyses)) {
		rasters_stack = stack(envVariables)
		names(rasters_stack) = envVariableNames1
		theRanges = ranges[h] # distance in meters
		gbm.x = gsub("\\.","-",names(rasters_stack))
		gbm.x = gbm.x[which(!gbm.x%in%envVariablesWithNotEnoughVariation)]
		gbm.y = "response"
		offset = NULL
		tree.complexity = 5 # "tc" = number of nodes in the trees
		learning.rate = 0.001 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(dataframes_list[[h]][[1]])[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "gaussian"
		n.trees = 10 # initial number of trees
		step.size = 2 # interval at which the predictive deviance is computed and logged
					  # (at each interval, the folds are successively used as test data set
					  # nd the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		keep.fold.vector = FALSE
		keep.fold.fit = FALSE
		showingFoldsPlot = FALSE
		brt_models_ccv = list() # classic cross-validations
		brt_models_scv = list() # spatial cross-validations
		for (i in 1:nberOfReplicates)
			{
				data = as.data.frame(dataframes_list[[h]][[i]]); n.trees = 10; learning.rate = 0.001; step.size = 2; fold.vector = NULL; worked = FALSE
				# BRT with classic (standard) cross-validation (CCV):
				pdf(file=paste0("All_BRT_models/2°_models/",subdirectory,"/BRT_model_",analyses[h],"_CCV_",i,".pdf"))
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								brt_models_ccv[[i]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
						if (length(brt_models_ccv) == i) worked = TRUE
					}
				dev.off()
				# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
				spdf = SpatialPointsDataFrame(data[c("x","y")], data[,c(1,4:dim(data)[2])], proj4string=crs(nullRaster)); worked = FALSE
				pdf(file=paste0("All_BRT_models/2°_models/",subdirectory,"/BRT_model_",analyses[h],"_SCV_",i,".pdf")); worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								myblocks = cv_spatial(spdf, column="response", k=n.folds, size=theRanges[1], selection="random")
								n.trees = 10; learning.rate = 0.001; step.size = 2; fold.vector = myblocks$folds_ids # "myblocks$foldID" not correct!
								brt_models_scv[[i]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
						if (length(brt_models_scv) == i) worked = TRUE
					}
				dev.off()
			}
		saveRDS(brt_models_ccv, paste0("All_BRT_models/2°_models/",subdirectory,"/BRT_models_",analyses[h],"_CCV.rds"))
		saveRDS(brt_models_scv, paste0("All_BRT_models/2°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds"))
	}}

# 5. Analysing the ecological niche models (RI, response curves, mapping)

analysis = 1 # analysis = 2
subdirectory = "BRT_models"

	# 5.1. Computation of the relative influence (RI) of each environmental factor

if (!file.exists(paste0(analysis,"°_RI_values.csv")))
	{
		relativeInfluences = matrix(0, nrow=length(envVariableNames1), ncol=3)
		row.names(relativeInfluences) = envVariableNames2; colnames(relativeInfluences) = paste0("RI_",analyses)
		for (h in 1:length(analyses))
			{
				brt_models_scv = readRDS(paste0("All_BRT_models/",analysis,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds"))
				for (i in 1:length(envVariableNames1))
					{
						allV = rep(NA, length(brt_models_scv))
						for (j in 1:length(brt_models_scv))
							{
								allV[j] = summary(brt_models_scv[[j]])[gsub("-","\\.",envVariableNames1)[i],"rel.inf"]
							}
						meanV = round(mean(allV),1); minV = round(min(allV),1); maxV = round(max(allV),1); quantiles = round(quantile(allV,c(0.025,0.975)),1)
						relativeInfluences[i,paste0("RI_",analyses[h])] = paste0(meanV," [",quantiles[1],"-",quantiles[2],"]")
					}
				for (i in 1:dim(relativeInfluences)[1])
					{
						txt = relativeInfluences[i,paste0("RI_",analyses[h])]; txts = rep(NA, 3)
						txts[1] = unlist(strsplit(txt," \\["))[1]
						temp = unlist(strsplit(txt," \\["))[2]
						txts[2] = unlist(strsplit(temp,"-"))[1]
						txts[3] = gsub("\\]","",unlist(strsplit(temp,"-"))[2])
						for (j in 1:length(txts))
							{
								if (!grepl("\\.",txts[j])) txts[j] = paste0(txts[j],".0")
							}
						relativeInfluences[i,paste0("RI_",analyses[h])] = paste0(txts[1]," [",txts[2],"-",txts[3],"]")
						if (grepl("NA",relativeInfluences[i,paste0("RI_",analyses[h])])) relativeInfluences[i,paste0("RI_",analyses[h])] = "-"
					}
			}
		write.table(relativeInfluences, paste0(analysis,"°_RI_values.csv"), quote=F, sep=",")
	}

	# 5.2. Generating the response curves for each environmental factor

for (h in 1:length(analyses))
	{
		dataframes_list = readRDS(paste0("Dataframes_",analysis,".rds"))
		envVariableValues = matrix(nrow=5, ncol=length(envVariableNames1))
		if (analysis == 1) 
			{
				envVariablesWithNotEnoughVariation = c("Perc_Bfor_ref_Grid_ETRS_LAEA","Perc_Cfor_ref_Grid_ETRS_LAEA")
			}
		if (analysis == 2) 
			{
				envVariablesWithNotEnoughVariation = c("Perc_Bfor_ref_Grid_ETRS_LAEA","Perc_Cfor_ref_Grid_ETRS_LAEA")
			}
		row.names(envVariableValues) = c("median","q1","q3","minV","maxV")
		colnames(envVariableValues) = envVariableNames1
		for (i in 1:length(envVariableNames1))
			{
				if (!envVariableNames1[i]%in%envVariablesWithNotEnoughVariation)
					{
						vS = c()
						for (j in 1:nberOfReplicates) vS = c(vS, dataframes_list[[h]][[j]][,gsub("-","\\.",envVariableNames1)[i]])
						minV = min(vS, na.rm=T); maxV = max(vS, na.rm=T); medianV = median(vS, na.rm=T)
						q1V = quantile(vS, probs=0.25, na.rm=T); q3V = quantile(vS, probs=0.75, na.rm=T)
						envVariableValues[,i] = cbind(medianV, q1V, q3V, minV, maxV)
					}
			}
		relativeInfluences = read.csv(paste0(analysis,"°_RI_values.csv"), head=T); quantiles = c("median"); options(scipen=10)
		brt_models_scv = readRDS(paste0("All_BRT_models/",analysis,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds"))
		if (savingPlots)
			{
				pdf(paste0(analysis,analyses[h],"_responses_NEW.pdf"), width=7, height=8.0) # dev.new(width=7, height=8.0)
				par(mfrow=c(6,5), oma=c(0.9,1.3,1,0.5), mar=c(2.75,1,0.5,1), lwd=0.2, bty="o", col="gray30", col.axis="gray30", fg="gray30")
				for (i in 1:length(envVariableNames1))
					{
						if (!envVariableNames1[i]%in%envVariablesWithNotEnoughVariation)
							{
								valuesInterval = 0.1; valuesInterval = (envVariableValues["maxV",i]-envVariableValues["minV",i])/100
								df = data.frame(matrix(nrow=length(seq(envVariableValues["minV",i],envVariableValues["maxV",i],valuesInterval)),ncol=length(envVariableNames1)))
								colnames(df) = gsub("-","\\.",envVariableNames1)
								for (j in 1:length(envVariableNames1))
									{
										valuesInterval = 0.1; valuesInterval = (envVariableValues["maxV",j]-envVariableValues["minV",j])/100
										if (i == j) df[,gsub("-","\\.",envVariableNames1)[j]] = seq(envVariableValues["minV",j],envVariableValues["maxV",j],valuesInterval)
										if (i != j) df[,gsub("-","\\.",envVariableNames1)[j]] = rep(envVariableValues[quantiles[1],j],dim(df)[1])
									}
								for (j in 1:length(brt_models_scv))
									{
										n.trees = brt_models_scv[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
										projection = predict.gbm(brt_models_scv[[j]], newdata=df, n.trees, type, single.tree)
										if ((i == 1)&(j == 1))
											{
												minY = min(projection); maxY = max(projection)
											}	else	{
												if (minY > min(projection)) minY = min(projection)
												if (maxY < max(projection)) maxY = max(projection)
											}
									}
							}
					}
				for (i in 1:length(envVariableNames1))
					{
						if (!envVariableNames1[i]%in%envVariablesWithNotEnoughVariation)
							{
								projections_list = list(); dfs = list()
								projections = list(); valuesInterval = 0.1; valuesInterval = (envVariableValues["maxV",i]-envVariableValues["minV",i])/100
								df = data.frame(matrix(nrow=length(seq(envVariableValues["minV",i],envVariableValues["maxV",i],valuesInterval)),ncol=length(envVariableNames1)))
								colnames(df) = gsub("-","\\.",envVariableNames1)
								for (j in 1:length(envVariableNames1))
									{
										valuesInterval = 0.1; valuesInterval = (envVariableValues["maxV",j]-envVariableValues["minV",j])/100
										if (i == j) df[,gsub("-","\\.",envVariableNames1)[j]] = seq(envVariableValues["minV",j],envVariableValues["maxV",j],valuesInterval)
										if (i != j) df[,gsub("-","\\.",envVariableNames1)[j]] = rep(envVariableValues[quantiles[1],j],dim(df)[1])
									}
								for (j in 1:length(brt_models_scv))
									{
										n.trees = brt_models_scv[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
										projection = predict.gbm(brt_models_scv[[j]], newdata=df, n.trees, type, single.tree)
										if (j == 1)
											{
												minX = min(df[,gsub("-","\\.",envVariableNames1)[i]]); maxX = max(df[,gsub("-","\\.",envVariableNames1)[i]])											
											}	else	{
												if (minX > min(df[,gsub("-","\\.",envVariableNames1)[i]])) minX = min(df[,gsub("-","\\.",envVariableNames1)[i]])
												if (maxX < max(df[,gsub("-","\\.",envVariableNames1)[i]])) maxX = max(df[,gsub("-","\\.",envVariableNames1)[i]])
											}
										projections[[j]] = projection
									}
								projections_list[[1]] = projections; dfs[[1]] = df; cols1 = list(); cols2 = list()
								cols1[[1]] = rgb(222,67,39,255,maxColorValue=255); cols2[[1]] = rgb(222,67,39,100,maxColorValue=255) # red
								cols1[[2]] = rgb(250,165,33,255,maxColorValue=255); cols2[[2]] = rgb(250,165,33,100,maxColorValue=255) # orange
								cols1[[3]] = rgb(70,118,187,255,maxColorValue=255); cols2[[3]] = rgb(70,118,187,100,maxColorValue=255) # blue
								for (j in 1:length(brt_models_scv))
									{
										if (j == 1)
											{
												plot(dfs[[1]][,gsub("-","\\.",envVariableNames1)[i]], projections_list[[1]][[j]], col=cols1[[1]], ann=F,
													 axes=F, lwd=0.2, type="l",, xlim=c(minX,maxX), ylim=c(minY,maxY+((maxY-minY)/3)))
											}	else	{
												lines(dfs[[1]][,gsub("-","\\.",envVariableNames1)[i]], projections_list[[1]][[j]], col=cols1[[1]], lwd=0.2)
											}
									}
								axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.040, col.axis="gray30", mgp=c(0,0.14,0))
								if ((analysis == 1)&(h == 1)) axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.040, col.axis="gray30", mgp=c(0,0.30,0), at=seq(0,1,0.2))
								if ((analysis == 1)&(h == 2)) axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.040, col.axis="gray30", mgp=c(0,0.30,0), at=seq(0,1,0.2))
								if ((analysis == 1)&(h == 3)) axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.040, col.axis="gray30", mgp=c(0,0.30,0), at=seq(0,1,0.1))
								# title(ylab="predicted values", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
								if (envVariableNames2[i] != "Density of close wild bird cases")
									{
										title(xlab=paste0(envVariableNames2[i]), cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
									}	else	{
										title(xlab=paste0(envVariableNames2[i],"   "), cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")
									}
								mtext(paste0("RI = ",relativeInfluences[i,paste0("RI_",analyses[h])]), side=3, line=-1.5, cex=0.55, col="gray30")
								box(lwd=0.2, col="gray30")
							}
					}
				dev.off()
			}
	}

	# 5.3. Projecting the risk map for the risk of local H5Nx circulation
	
for (h in 1:length(analyses))
	{
		brt_models_scv = readRDS(paste0("All_BRT_models/",analysis,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds"))
		rasters_stack = stack(envVariables); brt_projections = list(); saving_sd_rasters = FALSE
		df1 = as.data.frame(rasters_stack); df2 = df1[which(!is.na(rowMeans(df1))),]
		percentages = rep(NA, length(brt_models_scv)); CVs = rep(NA, length(brt_models_scv))
		for (i in 1:length(brt_models_scv))
			{
				n.trees = brt_models_scv[[i]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				projection = predict.gbm(brt_models_scv[[i]], df2, n.trees, type, single.tree)
				rast = rasters_stack[[1]]; rast[!is.na(rast[])] = projection; brt_projections[[i]] = rast
				maxV = max(rast[], na.rm=T); minV = min(rast[], na.rm=T); percentages[i] = ((maxV-minV)/minV)*100
				CVs[i] = (sd(rast[], na.rm=T)/mean(rast[], na.rm=T))*100
			}
		if (analyses[h] == "b")
			{
				meanV = round(mean(percentages),0); quantiles = round(quantile(percentages, c(0.025,0.975)),0)
				cat(meanV,"% (2.5th-97.5th percentiles range = [",quantiles[1],", ",quantiles[2],"])\n",sep="")
					# ecological suitability for the occurrence of at least one infected farm was on average 98% higher in 
					# the most suitable grid cell than in the least suitable one (95% range across SCV replicates: 7-262%).
			}
		cat("Mean CV = ",round(mean(CVs),0),"% [",round(quantile(CVs,0.025),0),", ",round(quantile(CVs,0.975),0),"]\n",sep="")
			# mean CV = 6% (2.5th-97.5th percentiles range = [1, 20]), 17% [1, 37] and 3% [0, 14] for analyses a, b, and c, respectively
		mean_projection = mean(stack(brt_projections)); stdv_projection = calc(stack(brt_projections), fun=sd)
		writeRaster(mean_projection, paste0("All_BRT_models/",analysis,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_mean.tif"), overwrite=T)
		writeRaster(stdv_projection, paste0("All_BRT_models/",analysis,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_stdv.tif"), overwrite=T)
	}

	# 5.4. Using the Boyce Index (BI) to assess the models on occurrence data collected in 2025 (no positive record in 2024 farm/hobby units)

new_data = read.csv("Data_210926.csv", head=T, sep=";"); new_data = new_data[which(new_data[,"final_date"]!=""),]
new_dates = decimal_date(ymd(gsub("/","-",gsub(" 00:00:00","",new_data[,"final_date"]))))
# new_data = new_data[which((new_dates>=2025)&(new_dates<2026)),]
# new_data = new_data[which((new_dates>=2026)&(new_dates<2027)),]
new_data = new_data[which((new_dates>=2025)&(new_dates<2027)),]
new_coordinates_lam = matrix(nrow=dim(new_data)[1], ncol=2); colnames(new_coordinates_lam) = c("longitude_Lambert","latitude_Lambert")
for (i in 1:dim(new_data)[1])
	{
		if ((!is.na(new_data[i,"longitude_all"]))&(new_data[i,"longitude_all"]!="Unknown"))
			{
				buffer = data.frame(lon=as.numeric(gsub(",",".",new_data[i,"longitude_all"])), lat=as.numeric(gsub(",",".",new_data[i,"latitude_all"])))
				coordinates(buffer) = c("lon", "lat"); crs(buffer) = CRS("+init=epsg:4326") # WGS84
				new_coordinates_lam[i,c("longitude_Lambert","latitude_Lambert")] = spTransform(buffer, crs(nullRaster))@coords
			}
	}
if (!"longitude_Lambert"%in%colnames(new_data)) new_data = cbind(new_data, new_coordinates_lam)
temp = new_data[which(is.na(raster::extract(nullRaster,new_data[,c("longitude_Lambert","latitude_Lambert")]))),]
new_data = new_data[which(!is.na(raster::extract(nullRaster,new_data[,c("longitude_Lambert","latitude_Lambert")]))),]
data = new_data[which(new_data[,"animal_type"]%in%c("Farm")),] # only computing the Boyce index for analyses "b" (farm outbreaks only)
positive_entries = which((data$Result_H5H7Other=="H5")&(data$Result_Pathotype=="HP")&(data$Result_AFLoccurence=="POSITIVE"))
positive_entries_coordinates = data[positive_entries,c("longitude_Lambert","latitude_Lambert")]
positive_entries_coordinates_text = paste0(positive_entries_coordinates[,1],"_",positive_entries_coordinates[,2])
positive_entries_coordinates_text_unique = unique(positive_entries_coordinates_text)
positive_entries_coordinates_unique = matrix(nrow=length(positive_entries_coordinates_text_unique), ncol=2)
for (i in 1:length(positive_entries_coordinates_text_unique))
	{
		positive_entries_coordinates_unique[i,] = unlist(strsplit(positive_entries_coordinates_text_unique[i],"_"))
	}
positive_entries_coordinates_unique = positive_entries_coordinates_unique[which(positive_entries_coordinates_unique[,1]!="NA"),]
positive_entries_coordinates_unique = as.data.frame(positive_entries_coordinates_unique)
positive_entries_coordinates_unique[,1] = as.numeric(positive_entries_coordinates_unique[,1])
positive_entries_coordinates_unique[,2] = as.numeric(positive_entries_coordinates_unique[,2])
positive_cells_unique = unique(cellFromXY(nullRaster, positive_entries_coordinates_unique))
presence_cells = xyFromCell(nullRaster, positive_cells_unique); # print(dim(presence_cells))
 
h = 2; brt_models_scv = readRDS(paste0("All_BRT_models/",analysis,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_SCV.rds"))
BIs = rep(NA,length(brt_models_scv)); bins_min = seq(0.0,0.9,0.1); bins_max = seq(0.1,1.0,0.1); brt_projections = list()
for (i in 1:length(brt_models_scv))
	{
		rasters_stack = stack(envVariables)
		df1 = as.data.frame(rasters_stack); df2 = df1[which(!is.na(rowMeans(df1))),]
		n.trees = brt_models_scv[[i]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
		projection = predict.gbm(brt_models_scv[[i]], df2, n.trees, type, single.tree)
		rast = rasters_stack[[1]]; rast[!is.na(rast[])] = projection; brt_projections[[i]] = rast
		vS_presence = raster::extract(rast, presence_cells); vS_all = rast[!is.na(rast[])]
		Fis = rep(NA, length(bins_min))
		for (j in 1:length(bins_min))
			{
				P = sum((vS_presence>=bins_min[j])&(vS_presence<bins_max[j]))/length(vS_presence)
				E = sum((vS_all>=bins_min[j])&(vS_all<bins_max[j]))/length(vS_all)
				Fis[j] = P/E
			}
		bins_mid = bins_min+0.05; bins_mid = bins_mid[which(!is.na(Fis))]; Fis = Fis[which(!is.na(Fis))]
		BIs[i] = cor(bins_mid, Fis, method="spearman")
	}
meanV = round(mean(BIs),2); quantiles = round(quantile(BIs,c(0.025,0.975)),2)
cat("\tmean BI = ",meanV,", 2.5th-97.5th percentiles range = [",quantiles[1],", ",quantiles[2],"]\n",sep="")
	# mean BI = 0.74, 2.5th-97.5th percentiles range = [0.04, 1.00] for 2025
	# mean BI = 0.91, 2.5th-97.5th percentiles range = [0.34, 1.00] for 2026
	# mean BI = 0.91, 2.5th-97.5th percentiles range = [0.87, 0.95] for 2025 and 2026 so far
	# plot(mean(stack(brt_projections))); points(presence_cells)
	# "We eventually used occurrence records from farms in 2025 and 2026 to date (June 29, 2026) to assess the capacity 
	# of our models to predict more recent cases, as no infected farms were reported in 2024. With a mean Boyce Index of
	# 0.91 (2.5th-97.5th percentiles range = [0.87, 0.95]), the models demonstrated a strong positive association between
	# predicted suitability and observed occurrence patterns, indicating that occurrence records were consistently 
	# concentrated in areas predicted to have higher ecological suitability."

# 6. Generating a composite figure with the different risk maps

	# Note: we should focus on Bernouilli analyses given that (i) the Gaussian distribution allows negative values
	#		and (ii) given that only three 1x1 km pixels have more than one (i.e. two) infected production units

gadm_0 = sf:::as_Spatial(st_transform(st_as_sf(shapefile("GADM_shapefiles/GADM_BEL_0.shp")),crs(nullRaster)))
gadm_2 = sf:::as_Spatial(st_transform(st_as_sf(shapefile("GADM_shapefiles/GADM_BEL_2.shp")),crs(nullRaster)))
gadm_2_sldt = as(gadm_2, "SpatialLines"); lines = list(); c = 0
for (i in 1:length(gadm_2_sldt)) # to only keep inland province borders
	{
		for (j in 1:length(gadm_2_sldt[i]@lines[[1]]@Lines))
			{
				coords1 = gadm_2_sldt[i]@lines[[1]]@Lines[[j]]@coords; indices = c()
				for (k in 1:dim(gadm_2_sldt[i]@lines[[1]]@Lines[[j]]@coords)[1])
					{
						present = FALSE
						for (l in 1:length(gadm_2_sldt))
							{
								if (i != l)
									{
										for (m in 1:length(gadm_2_sldt[l]@lines[[1]]@Lines))
											{
												coords2 = gadm_2_sldt[l]@lines[[1]]@Lines[[m]]@coords
												dists = sqrt(((coords1[k,1]-coords2[,1])^2)+((coords1[k,2]-coords2[,2])^2))
												if (min(dists) < 1000) present = TRUE
											}
									}
							}
						if (present) indices = c(indices, k)
					}
				if (gadm_2@data[i,"NAME_1"] == "Bruxelles")
					{
						c = c+1; line = Line(gadm_2_sldt[i]@lines[[1]]@Lines[[j]]@coords); lines[[c]] = Lines(list(line), ID=c)
					}	else	{
						line_coords = gadm_2_sldt[i]@lines[[1]]@Lines[[j]]@coords[indices,]
						if (dim(line_coords)[1] != 0)
							{
								lc = line_coords
								differences = sqrt(((lc[2:dim(lc)[1],1]-lc[1:dim(lc)[1]-1,1])^2)+((lc[2:dim(lc)[1],2]-lc[1:dim(lc)[1]-1,2])^2))
								if (max(differences) < 5000)
									{
										c = c+1; line = Line(line_coords); lines[[c]] = Lines(list(line), ID=c)
									}	else	{
										index = which(differences==max(differences))
										c = c+1; line = Line(line_coords[1:index,]); lines[[c]] = Lines(list(line), ID=c)
										c = c+1; line = Line(line_coords[(index+1):dim(line_coords)[1],]); lines[[c]] = Lines(list(line), ID=c)
									}
							}
					}
			}	
	}
sp_lines = SpatialLines(lines); gadm_2_sldt = sp_lines; crs(gadm_2_sldt) = crs(gadm_2)
colour_scale_1 = paste0(rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120)[10:110]),"CC")
colour_scale_2 = paste0(colorRampPalette(brewer.pal(9,"YlOrBr"))(120)[10:110],"CC")
farms_mean_list = list(); farms_stdv_list = list(); maxSDValue = 0.10
for (g in 1:1) # for (g in 1:2)
	{
		rast = raster(paste0("All_BRT_models/",g,"°_models/",subdirectory,"/BRT_models_",analyses[2],"_mean.tif"))
		index = which(envVariableNames1 == "Sciensano_farms_density"); farms_density = envVariables[[index]]
		rast[which(farms_density[]==0)] = NA; indices = which(!is.na(rast[]))
		farms = matrix(nrow=length(indices), ncol=3); colnames(farms) = c("x","y","estimate")
		farms[,c("x","y")] = xyFromCell(rast, indices); farms[,"estimate"] = raster::extract(rast, farms[,c("x","y")])
		farms = farms[order(farms[,"estimate"],decreasing=F),]; farms_mean_list[[g]] = farms
		rast = raster(paste0("All_BRT_models/",g,"°_models/",subdirectory,"/BRT_models_",analyses[2],"_stdv.tif"))
		index = which(envVariableNames1 == "Sciensano_farms_density"); farms_density = envVariables[[index]]
		rast[which(farms_density[]==0)] = NA; indices = which(!is.na(rast[]))
		farms = matrix(nrow=length(indices), ncol=3); colnames(farms) = c("x","y","estimate")
		farms[,c("x","y")] = xyFromCell(rast, indices); farms[,"estimate"] = raster::extract(rast, farms[,c("x","y")])
		farms = farms[order(farms[,"estimate"],decreasing=F),]; farms_stdv_list[[g]] = farms
	}
maxV = max(farms_mean_list[[1]][,"estimate"], na.rm=T); minV = min(farms_mean_list[[1]][,"estimate"], na.rm=T)
# round(((maxV-minV)/minV)*100,1) # areas the most at risk for the occurrence of a least one infected farm are 
		# 62% more at risk than the areas the least at risk for the occurrence of a least one infected farm

colour_scale_1_tab = cbind(seq(0.25,0.75,0.005), colour_scale_1)
colnames(colour_scale_1_tab) = c("ecological_suitability_value","HEX_colour_code")
write.csv(colour_scale_1_tab, "Colour_scale.csv", row.names=F, quote=F)

dev.new(width=8, height=4.5); plotting_farm_dots = FALSE; firstLayout = TRUE; firstLayout = FALSE
par(mfrow=c(2,3), oma=c(0,0,0,0.0), mar=c(0.0,0.0,0.0,0), lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
labels = rbind(cbind("A","C","B"),cbind("D","F","E")); plottingPresences = FALSE; plottingSD = TRUE; largeColourScale = TRUE
for (g in 1:2)
	{
		for (h in c(1,3,2))
			{
				if (plottingSD)
					{
						df2 = readRDS(paste0("Dataframes_1.rds"))[[2]][[1]]; farm_cases = df2[which(df2[,"response"]==1),c("x","y")]
						df3 = readRDS(paste0("Dataframes_1.rds"))[[3]][[1]]; hobby_cases = df3[which(df3[,"response"]==1),c("x","y")]
						rast = raster(paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_models_",analyses[h],"_mean.tif"))
					}	else	{
						df2 = readRDS(paste0("Dataframes_",g,".rds"))[[2]][[1]]; farm_cases = df2[which(df2[,"response"]==1),c("x","y")]
						df3 = readRDS(paste0("Dataframes_",g,".rds"))[[3]][[1]]; hobby_cases = df3[which(df3[,"response"]==1),c("x","y")]
						rast = raster(paste0("All_BRT_models/",g,"°_models/",subdirectory,"/BRT_models_",analyses[h],"_mean.tif"))						
					}
				colour_scale_mod = colour_scale_1
				if (largeColourScale == TRUE)
					{
						index1 = round((((min(rast[], na.rm=T)-0.25)/(0.75-0.25))*100)+1)
						index2 = round((((max(rast[], na.rm=T)-0.25)/(0.75-0.25))*100)+1)
						colour_scale_mod = colour_scale_1[index1:index2]
					}
				if ((g == 2)&(plottingSD))
					{
						rast = raster(paste0("All_BRT_models/1°_models/",subdirectory,"/BRT_models_",analyses[h],"_stdv.tif"))						
						index1 = round((((min(rast[], na.rm=T)-0)/(maxSDValue-0))*100)+1)
						index2 = round((((max(rast[], na.rm=T)-0)/(maxSDValue-0))*100)+1)
						colour_scale_mod = colour_scale_2[index1:index2]
					}
				if (analyses[h] == "b")
					{
						farms_density = envVariables[[which(envVariableNames1 == "Sciensano_farms_density")]]; rast[which(farms_density[]==0)] = NA
					}
				plot(contour, col="gray95", border=NA, ann=F)
				if (analyses[h] != "b") { plot(rast, col="gray50", add=T, legend=F); plot(rast, col=colour_scale_mod, add=T, legend=F) }
				if (analyses[h] != "b") plot(gadm_2_sldt, col="gray40", lwd=0.4, add=T)
				if (analyses[h] == "b") { plot(rast, col="gray50", add=T, legend=F); plot(rast, col=colour_scale_mod, add=T, legend=F) }
				if (analyses[h] == "b") plot(gadm_2_sldt, col="gray40", lwd=0.4, add=T)
				plot(contour, col=NA, border="gray40", lwd=0.5, add=T)
				if ((h == 1)&(plottingPresences))
					{
						points(farm_cases, cex=0.45, col="black", lwd=0.2, pch=16)
						points(hobby_cases, cex=0.5, col="gray30", lwd=0.2, pch=17)
					}
				if ((h == 2)&(plotting_farm_dots))
					{
						if (!((g == 2)&(plottingSD)))
							{
								vS = farms_mean_list[[g]][,"estimate"]; cols = colour_scale_mod[(((vS-min(vS))/(max(vS)-min(vS)))*100)+1]
								if (largeColourScale == TRUE) cols = colour_scale_1[(((vS-0.25)/(0.75-0.25))*100)+1]
								for (i in 1:dim(farms_mean_list[[g]])[1])
									{
										points(farms_mean_list[[g]][i,"x"], farms_mean_list[[g]][i,"y"], cex=0.45, col="gray50", pch=16)
										points(farms_mean_list[[g]][i,"x"], farms_mean_list[[g]][i,"y"], cex=0.45, col=cols[i], pch=16)
										points(farms_mean_list[[g]][i,"x"], farms_mean_list[[g]][i,"y"], cex=0.45, col="gray30", lwd=0.1, pch=1)
									}
							}	else	{
								vS = farms_stdv_list[[g]][,"estimate"]; cols = colour_scale_2[(((vS-0)/(maxSDValue-0))*100)+1]
								for (i in 1:dim(farms_stdv_list[[g]])[1])
									{
										points(farms_stdv_list[[g]][i,"x"], farms_stdv_list[[g]][i,"y"], cex=0.45, col="gray50", pch=16)
										points(farms_stdv_list[[g]][i,"x"], farms_stdv_list[[g]][i,"y"], cex=0.45, col=cols[i], pch=16)
										points(farms_stdv_list[[g]][i,"x"], farms_stdv_list[[g]][i,"y"], cex=0.45, col="gray30", lwd=0.1, pch=1)
									}								
							}
					}
				buffer = rast; colour_scale_legend = colour_scale_1
				if ((g == 2)&(plottingSD)) colour_scale_legend = colour_scale_2
				par(lwd=0.1, bty="o", col="gray30", col.axis="gray30", fg=NA)
				if (firstLayout)
					{
						plot(rast, legend.only=T, add=T, col="gray50", legend.width=0.5, legend.shrink=0.3, col.tick=NA,
							 smallplot=c(0.870,0.885,0.055,0.270), legend.args=list(text="", cex=0.7, line=0.0, col=NA, col.axis=NA), horizontal=F,
							 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0, tck=-1.0, col=NA, col.tick=NA, col.axis=NA, mgp=c(0,0.5,0)))
						par(lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
						if (largeColourScale == TRUE) buffer = raster(as.matrix(cbind(0.25,0.75)))
						if ((g == 2)&(plottingSD)) buffer = raster(as.matrix(cbind(0,maxSDValue)))
						plot(buffer, legend.only=T, add=T, col=colour_scale_legend, legend.width=0.5, legend.shrink=0.3,
							 smallplot=c(0.870,0.885,0.055,0.270), legend.args=list(text="", cex=0.7, line=0.0, col="gray30"), horizontal=F,
					 		 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0.3, tck=-1.0, col.tick="gray30", col.axis="gray30", mgp=c(0,0.5,0)))
						mtext(labels[g,h], side=3, line=-1.7, at=3800000, cex=0.8, col="gray30", adj=0)
						mtext("Risk mapping for H5Nx", side=3, line=-14.4, at=3827000, cex=0.55, col="gray30", adj=0)
						if (h == 1) mtext("occurrence in hobby/farms", side=3, line=-15.3, at=3827000, cex=0.55, col="gray30", adj=0)
						if (h == 2) mtext("occurrence in farms", side=3, line=-15.3, at=3827000, cex=0.55, col="gray30", adj=0)
						if (h == 3) mtext("occurrence in hobby", side=3, line=-15.3, at=3827000, cex=0.55, col="gray30", adj=0)
						if (g == 1) mtext("(variable: min. 1 infected unit)", side=3, line=-16.2, at=3827000, cex=0.55, col="gray30", adj=0)
						if ((g == 2)&(!plottingSD)) mtext("(variable: infected units/pixel)", side=3, line=-16.2, at=3827000, cex=0.55, col="gray30", adj=0)
						if ((g == 2)&(plottingSD)) mtext("(standard deviation)", side=3, line=-16.2, at=3827000, cex=0.55, col="gray30", adj=0)
					}	else	{
						plot(rast, legend.only=T, add=T, col="gray50", legend.width=0.5, legend.shrink=0.3, col.tick=NA,
							 smallplot=c(0.043,0.058,0.26,0.57), legend.args=list(text="", cex=0.7, line=0.0, col=NA, col.axis=NA), horizontal=F,
							 axis.args=list(cex.axis=0.65, lwd=0, lwd.tick=0, tck=-1.0, col=NA, col.tick=NA, col.axis=NA, mgp=c(0,0.5,0)))
						par(lwd=0.3, bty="o", col="gray30", col.axis="gray30", fg="gray30")
						if (largeColourScale == TRUE) buffer = raster(as.matrix(cbind(0.25,0.75)))
						if ((g == 2)&(plottingSD)) buffer = raster(as.matrix(cbind(0,maxSDValue)))
						plot(buffer, legend.only=T, add=T, col=colour_scale_legend, legend.width=0.5, legend.shrink=0.3,
							 smallplot=c(0.043,0.058,0.26,0.57), legend.args=list(text="", cex=0.7, line=0.0, col="gray30"), horizontal=F,
					 		 axis.args=list(cex.axis=0.85, lwd=0, lwd.tick=0.3, tck=-1.0, col.tick="gray30", col.axis="gray30", mgp=c(0,0.5,0)))
						mtext(labels[g,h], side=3, line=-1.7, at=3800000, cex=0.8, col="gray30", adj=0)
						mtext("Risk mapping for H5Nx", side=3, line=-14.4, at=3800000, cex=0.65, col="gray30", adj=0)
						if (h == 1) mtext("occurrence in hobby/farms", side=3, line=-15.4, at=3800000, cex=0.65, col="gray30", adj=0)
						if (h == 2) mtext("occurrence in farms", side=3, line=-15.4, at=3800000, cex=0.65, col="gray30", adj=0)
						if (h == 3) mtext("occurrence in hobby", side=3, line=-15.4, at=3800000, cex=0.65, col="gray30", adj=0)
						if (g == 1) mtext("(variable: min. 1 infected unit)", side=3, line=-16.4, at=3800000, cex=0.65, col="gray30", adj=0)
						if ((g == 2)&(!plottingSD)) mtext("(variable: infected units/pixel)", side=3, line=-16.4, at=3800000, cex=0.55, col="gray30", adj=0)
						if ((g == 2)&(plottingSD)) mtext("(standard deviation)", side=3, line=-16.4, at=3800000, cex=0.65, col="gray30", adj=0)						
					}
			}
	}
# system(paste0("magick -units PixelsPerInch -density 1000 Risk_maps_1.pdf -background white -alpha remove -flatten Risk_maps_1.png"))
# system(paste0("magick -units PixelsPerInch -density 1000 Risk_maps_2.pdf -background white -alpha remove -flatten Risk_maps_2.png"))

