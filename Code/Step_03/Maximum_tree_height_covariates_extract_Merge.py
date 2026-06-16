from os import listdir
from os.path import isfile, join

# fout=open("CovariatesExtractionFolder/20201012_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel.csv","a")
fout=open("CovariatesExtractionFolder/Tree_height_covariates_extraction_merged_20260407_filtered.csv","a")

sampled_data = listdir("CovariatesExtractionFolder/SubExtractedTables")
sampled_data = list(filter(lambda f: f.endswith('.csv'), sampled_data))

# first file:
for line in open("CovariatesExtractionFolder/SubExtractedTables/"+sampled_data[0]):
    fout.write(line)
# now the rest:
for num in range(1,len(sampled_data)):
	try:
		f = open("CovariatesExtractionFolder/SubExtractedTables/"+sampled_data[num])
		f.__next__()
		for line in f:
			fout.write(line)
		f.close()
	except IOError:
	    pass

fout.close()
