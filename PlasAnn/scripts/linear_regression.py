"""
LINEAR_REGRESSION.PY

DESCRIPTION: Runs regressions on using files created from LR_FILE.PY
"""
import pandas as pd
import numpy as np
import itertools
import scipy.stats as scipystats
import statsmodels.api as sm
from statsmodels.graphics.regressionplots import *
import matplotlib.pyplot as plt
import xlsxwriter as xl
from glob import glob
import os
import shutil

#Create folder
os.system("mkdir ./../output/linear-reg-out")

ch_path = "./../lr_files/"
folder = glob("%s/*.csv" %ch_path)

if not os.path.exists('./../lr_files/figures/'):
	os.makedirs('./../lr_files/figures/')

if not os.path.exists('./../output/linear-reg-out/top10/'):
	os.makedirs('./../output/linear-reg-out/top10/')
for file in folder: 

	regress = file[11:-4]

	#Read in CSV File
	df = pd.read_csv(file)
	func_df = pd.DataFrame()
	if "cog" in regress: 
		func_df = pd.read_csv("./../lr_files/functions/cog_functions.csv")
	elif "go" in regress: 
		func_df = pd.read_csv("./../lr_files/functions/go_functions.csv")

	results_df = pd.DataFrame(columns = ["Reg_Form",  "R^2", "Adjusted R^2", "F-statistic", "P-Val of F-Statistic", "Slope", "Intercept", "Function"])

	#Open results.txt File
	append=0
	if os.path.exists("./../output/linear-reg-out/"+regress +"_results.txt"):
		append=1
	results=open("./../output/linear-reg-out/"+regress +"_results.txt", "a")
	if append==1:
		results.write("\n\n***** APPENDING NEW RUN *****")

	y = df["AC_avg"]

	# Removing characters from columns that cause problems in OLS 
	rem_char = ['(', ')', '.', '-', '/', '\'']
	for column in df.columns: 
		new_col = column
		if "|" not in column: 
			for i in rem_char: 
				new_col = new_col.replace(i, '_')
			df = df.rename({column : new_col}, axis = 1)

	plasmid_features = df.iloc[0:df.shape[0],6:]
	print(plasmid_features.head())

	res_sum = []
	func = ""

	for (columnName, columnData) in plasmid_features.iteritems():
		X = plasmid_features[columnName]
		
		if  "|" not in columnName and np.count_nonzero(X) > 2: 
			print("\n===" + str(columnName) + "===")
			max_X = X.max()

			if regress != "cog_func" or regress != "gene_absence_presence":
				if columnName in func_df.columns: 
					func = str(func_df.at[0,columnName])
					results.write("\n" + str(columnName) + ": " + func + "\n")
				else: 
					results.write("\n" + str(columnName) + "\n")
					func = ""
			x = sm.add_constant(X)
			reg = sm.OLS(y, x).fit()
			print(reg.summary())

			form = str(columnName).replace('/', '')+" ~ AC_avg"
			p = reg.params
			const = p[0]
			slope = p[1]

			test = np.linspace(-1,max_X+1)
			ax = df.plot(x = columnName, y = "AC_avg", kind = "scatter")
			ax.plot(test, const + slope*test)
			ax.set_xlabel(columnName +" | "+func, fontsize = 18)
			ax.set_ylabel("AC_avg", fontsize = 18)
			ax.set_ylim([0,2])
			ax.set_xlim([-1,max_X+1])

			fig = ax.figure

			fig.savefig('./../lr_files/figures/'+columnName.replace('/', '')+'.png', dpi=fig.dpi)

			res_sum = [form, reg.rsquared, reg.rsquared_adj, reg.fvalue, reg.f_pvalue, const, slope, func]

			results_df.loc[len(results_df)] = res_sum
			
			results.write(str(reg.summary()))

	results_df.sort_values(by = ["R^2", "P-Val of F-Statistic"], inplace = True, ascending = [False, True])
	results_df = results_df.round({ "R^2":2,  "Adjusted R^2":2, "F-statistic":2, "P-Val of F-Statistic":2, "Slope":2, "Intercept":2})
	results_df.to_csv("./../output/linear-reg-out/"+regress.replace('/','') +"_results.csv")
#----- 
	count = 0
	df_top10 = results_df.iloc[:10]
	for index, row in df_top10.iterrows(): 


		ind = row["Reg_Form"].index("~")
		col = str(row["Reg_Form"])[:ind-1]

		file = './../lr_files/figures/'+col.replace('/', '')+'.png'
		shutil.move(file, './../output/linear-reg-out/top10/'+col.replace('/', '')+'.png')

#----- 
shutil.rmtree('./../lr_files/figures/')
