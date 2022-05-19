# Code Used from https://towardsdatascience.com/how-to-build-a-regression-model-in-python-9a10685c7f09
# AND https://songhuiming.github.io/pages/2017/01/21/linear-regression-in-python-chapter-3-regression-with-categorical-predictors/
#
#
# Creates Linear Regression Models for given plasmid data + acquisition cost data
#

import pandas as pd
import numpy as np
import itertools
from itertools import chain, combinations
import statsmodels.formula.api as smf
import scipy.stats as scipystats
import statsmodels.api as sm
import statsmodels.stats.stattools as stools
import statsmodels.stats as stats
from statsmodels.graphics.regressionplots import *
import matplotlib.pyplot as plt
import xlsxwriter as xl
import copy
import math
import time
import os

#Create folder
os.system("mkdir ./../output/linear-reg-out")

#Read in CSV File
df = pd.read_csv("./../output/PlasmidRegression.csv")
print(df)

#Open File
append=0
if os.path.exists("./../output/linear-reg-out/results.txt"):
	append=1
results=open("./../output/linear-reg-out/results.txt", "a")
if append==1:
	results.write("\n\n***** APPENDING NEW RUN *****")

#Remove unecessary columns
df=df.drop(columns=['Number', 'Replication-Genes','Conjugation-Genes', 'Metabolism-Genes', 'Stress-Genes', 'Toxin-Genes', 'Antitoxin-Genes', 'AC_std'])

#Seperate gene presence data into its own df
data_index=df.columns.get_loc("Antibiotic")
gene_pres=df.iloc[0:df.shape[0]+1,2:data_index]
print(gene_pres)

#Create xlsx file for
labels=["Reg_Form", "R^2", "Adjusted R^2", "F-statistic", "P-Val of F-Statistic"]
workbook=xl.Workbook('./../output/linear-reg-out/lr_results.xlsx')
worksheet=workbook.add_worksheet()
worksheet.write_row('A1', labels)

#---Linear Regression Start---
#Based on length
print("\n\n--Length, no covariates, intercept==0--")
reg=smf.ols(formula="normalized_AC ~ 0 + Length", data=df).fit()
print(reg.summary())
results.write("---Length, no covariates, intercept=0---\n")
results.write(str(reg.summary()))
worksheet.write_row('A1', ['normalized_AC ~ 0 + Length', str(reg.rsquared), str(reg.rsquared_adj),str(reg.fvalue), str(reg.f_pvalue)])

#Based on Antibiotic Gene Numbers
print("\n\n--Antibiotic Genes, no covariates, intercept==0--")
reg=smf.ols(formula="normalized_AC ~ 0 + Antibiotic-Genes", data=df).fit()
print(reg.summary())
results.write("---Antibiotic Genes, no covariates, intercept=0---\n")
results.write(str(reg.summary()))
worksheet.write_row('A2', ['normalized_AC ~ 0 + Antibiotic-Genes', str(reg.rsquared), str(reg.rsquared_adj),str(reg.fvalue), str(reg.f_pvalue)])

#Based on INC group, no covariates, with intercept == 0
print("\n\n--INC Group, no covariates, intercept==0--")
reg=smf.ols(formula="normalized_AC ~ 0 + INC", data=df).fit()
print(reg.summary())
results.write("---INC Group, no covariates, intercept=0---\n")
results.write(str(reg.summary()))
worksheet.write_row('A3', ['normalized_AC ~ 0 + INC', str(reg.rsquared), str(reg.rsquared_adj),str(reg.fvalue), str(reg.f_pvalue)])

print("\n\n---Based on Antibiotic, no covariates---")
reg=smf.ols(formula="normalized_AC ~ 0 + Antibiotic", data=df).fit()
print(reg.summary())
results.write("---Based on Antibiotic, no covariates---\n")
results.write(str(reg.summary()))
worksheet.write_row('A4', ['normalized_AC ~ 0 + Antibiotic', str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])

print("\n\n---Use Carb as reference---")
reg=smf.ols(formula="normalized_AC ~ Antibiotic", data=df).fit()
results.write("---Use Carb as reference---\n")
results.write(str(reg.summary()))
print(str(reg.summary()))
worksheet.write_row('A5', ['normalized_AC ~ Antibiotic', str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])

print("\n\n---Use Antibiotic with INC covariates---")
reg=smf.ols(formula="normalized_AC ~ 0 + Antibiotic + INC", data=df).fit()
results.write("---Use Antibiotic with INC covariates---\n")
results.write(str(reg.summary()))
print(str(reg.summary()))
worksheet.write_row('A6', ['normalized_AC ~ 0 + Antibiotic + INC', str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])

#Holds full formula for ALL VARIABLES LR analysis
full_form="normalized_AC ~ 0 + Antibiotic + INC"

print("\n\n---Gene Presence Data Added---")
results.write("---Gene Presence Data Added---\n")
xl_index=6
for (columnName, columnData) in gene_pres.iteritems():
	print("\n===" + str(columnName) + "===")
	results.write("\n===" + str(columnName) + "===")
	full_form=full_form + ' + Q(' + '"' + columnName + '"' + ')'

	print("\n---Gene Presence Only---")
	results.write("\n---Gene Presence Only---\n")
	form='normalized_AC ~ 0 + Q(' + '"' + columnName + '"' + ')'
	reg=smf.ols(formula=form, data=df).fit()
	print(reg.summary())
	results.write(str(reg.summary()))
	worksheet.write_row('A' + str(xl_index), [form, str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])
	xl_index=xl_index+1

	print("\n---Gene Presence and INC Group---")
	results.write("\n---Gene Presence and INC Group---\n")
	form='normalized_AC ~ 0 + INC + Q(' + '"' + columnName + '"' + ')'
	reg=smf.ols(formula=form, data=df).fit()
	print(reg.summary())
	results.write(str(reg.summary()))
	worksheet.write_row('A' + str(xl_index), [form, str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])
	xl_index=xl_index+1

	print("\n---Gene Presence and Antibiotic---")
	results.write("\n---Gene Presence and Antibiotic---\n")
	form='normalized_AC ~ 0 + Antibiotic + Q(' + '"' + columnName + '"' + ')'
	reg=smf.ols(formula=form, data=df).fit()
	print(reg.summary())
	results.write(str(reg.summary()))
	worksheet.write_row('A' + str(xl_index), [form, str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])
	xl_index=xl_index+1

	print("\n---Gene Presence, INC, and Antibiotic---")
	results.write("\n---Gene Presence, INC, and Antibiotic---\n")
	form='normalized_AC ~ 0 + Antibiotic + INC + Q(' + '"' + columnName + '"' + ')'
	reg=smf.ols(formula=form, data=df).fit()
	print(reg.summary())
	results.write(str(reg.summary()))
	worksheet.write_row('A' + str(xl_index), [form, str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])
	xl_index=xl_index+1

print("\n\n---ALL VARIABLES---")
reg=smf.ols(formula=full_form, data=df).fit()
print(str(reg.summary()))
results.write("\n---ALL VARIABLES---\n")
results.write(str(reg.summary()))
worksheet.write_row('A' + str(xl_index), ['ALL VARIABLES', str(reg.rsquared), str(reg.rsquared_adj), str(reg.fvalue), str(reg.f_pvalue)])

workbook.close()
