all: submission

submission:
	# Some directories we need
	mkdir training missing_trainingchunk_predictions imputedDataWithFactorLoadings sample_code
	# Compute month, hour predictions
	python/median_hour_month.py > missing_trainingchunk_predictions/median_3month.csv
	# Run sample code provided for baseline
	Rscript r/sample_code.r
	# Impute data / reduce factors with SVD...
	# split_chunks - create training and test sets for imputed data
	# prep_predictions 
	# blending

clean: 
	rm -rf training missing_trainingchunk_predictions imputedDataWithFactorLoadings sample_code