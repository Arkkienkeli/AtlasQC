import csv
import os
import random
from  math import atan2
from math import atan
from math import pi
from math import degrees
from matplotlib.lines import Line2D
import numpy as np
from sklearn import linear_model
from scipy.spatial.distance import euclidean
import argparse
import subprocess
import sqlite3
from itertools import groupby
import collections
import time
import cProfile
import plotly.plotly as py
import plotly.graph_objs as go
import plotly

class Forel(object):

	def __init__(self, angle, X, key, refgenotypes, chromosome, gender, real_genotypes, X_Raw):
		self.X_dict = X
		self.X_dict_Raw = X_Raw
		self.chromosome = chromosome
		self.key = key
		self.result = []
		self.resultR = []
		self.resultRaw = []
		self.resultG = []
		self.resultMT = []
		self.resultXH = []
		self.resultYH = []
		self.resultYF = []
		self.genotypes = {}
		self.refgenotypes = refgenotypes
		self.real_genotypes = real_genotypes
		self.flag = 0
		self.gender = gender
		self.R = angle
		self.HW = True
		self.resultCluster = []
		self.main()


	def main(self):

		mas = []
		classes = []
		centers = []
		if self.R == None:
			self.R = 10.0

		X = []
		X_Raw = []
		labels = []

		for i in range(len(self.X_dict)):
			X_Raw.append(list(self.X_dict_Raw[i].items())[0][1])
			X.append(list(self.X_dict[i].items())[0][1])
			labels.append(list(self.X_dict[i].items())[0][0])


		if self.gender == None:
			self.gender = {}
			for l in labels:
				self.gender[l] = 'U'

		xmax_party = 0
		for i in range(len(X)):
			if (X[i][0] > 0.200 or X[i][1] > 0.200) and (X_Raw[i][0] > 500 or X_Raw[0][1] > 500):
				mas.append(X[i])
			elif (X[i][0] > 0.200 or X[i][1] > 0.200) and (X_Raw[i][0] < 500 and X_Raw[0][1] < 500):
				self.result.append(labels[i])
				self.resultRaw.append(labels[i])
			elif (X[i][0] < 0.200 and X[i][1] < 0.200) and (X_Raw[i][0] < 500 and X_Raw[0][1] < 500):
				self.result.append(labels[i])
				self.resultR.append(labels[i])
				self.resultRaw.append(labels[i])
			else:
				self.result.append(labels[i])
				self.resultR.append(labels[i])
			if X[i][0] > xmax_party:
				xmax_party = X[i][0]


		k = 0
		while(mas != []):
			current = degrees(np.arctan2(mas[0][1], mas[0][0]))
			neighbours = self.get_neighbours(mas, current, self.R)
			center = self.get_center(neighbours)
			while(center != current and center != current):
				current = center
				neighbours = self.get_neighbours(mas, current, self.R)
				center = self.get_center(neighbours)
				if k > 3:
					break
				k = k + 1
			neighboursNew = self.get_neighbours(mas, current, self.R)
			classes.append(neighboursNew)
			centers.append(center)
			mas = self.delete_objs(mas, neighbours)

		degrees_ = []
		lrr = []

		failed_clusters = []
		cental_line_angles = []
		central_angles = []


		xy = []
		for i in range(len(classes)):
			xmax = 0.000000000001
			ymax = 0
			for j in range(len(classes[i])):
				if classes[i][j][0] > xmax:
					xmax = classes[i][j][0]
				if classes[i][j][1] > ymax:
					ymax = classes[i][j][1]

			y_std = np.std(np.array(list(zip(*classes[i]))[1]))
			x_std = np.std(np.array(list(zip(*classes[i]))[0]))

			x_range = max(np.array(list(zip(*classes[i]))[0])) - min(np.array(list(zip(*classes[i]))[0]))
			y_range = max(np.array(list(zip(*classes[i]))[1])) - min(np.array(list(zip(*classes[i]))[1]))
			if x_std == 0.000000:
				x_std = 0.00000001
			if (ymax / xmax) > 2:
				lr = linear_model.LinearRegression().fit(np.array(list(zip(*classes[i]))[1]).reshape(-1, 1), np.array(list(zip(*classes[i]))[0]))
				lr_result = lr.predict(np.array(list(zip(*classes[i]))[1]).reshape(-1, 1))
				if (len(classes[i]) == 1):
					angle = centers[i]
					angle_center = np.rad2deg(np.arctan2(np.array(list(zip(*classes[i]))[1])[0], np.array(list(zip(*classes[i]))[0])[0]))
				else:
					angle = centers[i]
					angle_center = np.rad2deg(np.arctan2( ( abs(np.array(list(zip(*classes[i]))[1])[-1] - np.array(list(zip(*classes[i]))[1])[0])), abs((lr_result[-1] - lr_result[0]))))

				cental_line_angles.append(angle_center)
				central_angles.append(angle)
				lrr.append(lr)
				xy.append(1)
			else:
				lr = linear_model.LinearRegression().fit(np.array(list(zip(*classes[i]))[0]).reshape(-1, 1), np.array(list(zip(*classes[i]))[1]))
				lr_result = lr.predict(np.array(list(zip(*classes[i]))[0]).reshape(-1, 1))
				if (len(classes[i]) == 1):
					angle = centers[i]
					angle_center = np.rad2deg(np.arctan2(np.array(list(zip(*classes[i]))[1])[0], np.array(list(zip(*classes[i]))[0])[0]))
				else:
					angle = centers[i]
					angle_center = np.rad2deg(np.arctan2( abs(lr_result[-1] - lr_result[0]), abs(np.array(list(zip(*classes[i]))[0])[-1] - np.array(list(zip(*classes[i]))[0])[0])))


				cental_line_angles.append(angle_center)
				central_angles.append(angle)
				lrr.append(lr)
				xy.append(0)


		#some shit code to mark eliminated points by R
		points = []
		for i in range(len(classes)):
			for j in range(len(classes[i])):
				points.append(classes[i][j])

		#detecet heterozgosity
		#AA 0 AB 1 BB 2
		heterozgosity_cluster = []
		for i in range(len(classes)):
			if xy[i] == 1:
				heterozgosity_cluster.append(0)
			else:
				xmax, ymax = 0, 0.0000001
				for j in range(len(classes[i])):
					if classes[i][j][0] > xmax:
						xmax = classes[i][j][0]
					if classes[i][j][1] > ymax:
						ymax = classes[i][j][1]
				if xmax / ymax > 2:
					heterozgosity_cluster.append(2)
				else:
					heterozgosity_cluster.append(1)

		for i in range(len(classes)):
			if (heterozgosity_cluster[i] == 0 and central_angles[i] > 75):
				failed_clusters.append(0)
			elif (heterozgosity_cluster[i] == 1 and central_angles[i] > 30 and central_angles[i] < 60):
				failed_clusters.append(0)
			elif (heterozgosity_cluster[i] == 2 and central_angles[i] > 0 and central_angles[i] < 15):
				failed_clusters.append(0)
			else:
				failed_clusters.append(1)


		#get biggest cluster for each het-y state (but not AB)
		states = [0,2]
		clusters_count_arr = []
		for s in range(len(states)):
			clusters_count = 0
			for i in range(len(classes)):
				if failed_clusters[i] == 0 and heterozgosity_cluster[i] == states[s]:
					clusters_count = clusters_count + 1
				clusters_count_arr.append(clusters_count)
			if clusters_count > 1:
				maxCl = 0
				for i in range(len(classes)):
					if failed_clusters[i] == 0 and heterozgosity_cluster[i] == states[s] and len(classes[i]) > maxCl:
						maxCl =  len(classes[i])

				maxClclusters = []
				for i in range(len(classes)):
					if failed_clusters[i] == 0 and heterozgosity_cluster[i] == states[s] and len(classes[i]) < maxCl:
						failed_clusters[i] = 1
					elif failed_clusters[i] == 0 and heterozgosity_cluster[i] == states[s] and len(classes[i]) == maxCl:
						maxClclusters.append(i)

				if len(maxClclusters) > 1:
					mainCluster = maxClclusters[0]
					for i in range(1, len(maxClclusters)):
						for p in classes[maxClclusters[i]]:
							classes[mainCluster].append(p)
							for i in range(1, len(maxClclusters)):
								classes[maxClclusters[i]] = []
								failed_clusters[maxClclusters[i]] = 1

		
		ymax = 1.0
		xmax = 1.0
		for i in range(len(X)):
			if X[i][0] > xmax:
				xmax = X[i][0]

			if X[i][1] > ymax:
				ymax = X[i][1]

		for i in range(len(classes)):
			if failed_clusters[i] == 0 and heterozgosity_cluster[i] == 0:
				distance = 0.03 * xmax
				cluster_max = max(list(zip(*classes[i]))[0])
				for j in range(len(classes)):
					if failed_clusters[j] == 1 and heterozgosity_cluster[j] == 0:
						for k in range(len(classes[j])):
							if classes[j][k][0] <= cluster_max + distance:
								classes[i].append(classes[j][k])
						classes[j] = [x for x in classes[j] if x[0] > cluster_max + distance]


		for i in range(len(classes)):
			if failed_clusters[i] == 0 and heterozgosity_cluster[i] == 2:
				distance = 0.03 * ymax
				cluster_max = max(list(zip(*classes[i]))[1])
				for j in range(len(classes)):
					if failed_clusters[j] == 1 and heterozgosity_cluster[j] == 2:
						for k in range(len(classes[j])):
							if classes[j][k][1] <= cluster_max + distance:
								classes[i].append(classes[j][k])
						classes[j] = [x for x in classes[j] if x[1] > cluster_max + distance]

		for i in range(len(classes)):
			if failed_clusters[i] == 1:
				for j in range(len(X)):
					if X[j] in classes[i]:
						self.resultCluster.append(labels[j])
						self.result.append(labels[j])
		
		cl = []
		ymax = 1.0
		xmax = 1.0

		for i in range(len(X)):
			if X[i][0] > xmax:
				xmax = X[i][0]

			if X[i][1] > ymax:
				ymax = X[i][1]

		self.result = list(set(self.result))

		notfailed_hetero_clusters = [heterozgosity_cluster[i] for i in range(len(classes)) if failed_clusters[i] == 0]
		if 0 in notfailed_hetero_clusters and 2 in notfailed_hetero_clusters and 1 not in notfailed_hetero_clusters:
			self.flag = self.flag + 2


		for i in range(len(classes)):
			if (heterozgosity_cluster[i] == 0 or heterozgosity_cluster[i] == 2) and failed_clusters[i] == 0:
				for j in range(len(X)):
					if X[j] in classes[i]:
						self.genotypes[labels[j]] = 'AA'
			elif heterozgosity_cluster[i] == 1  and failed_clusters[i] == 0:
				for j in range(len(X)):
					if X[j] in classes[i]:
						self.genotypes[labels[j]] = 'AB'


		for res in self.genotypes.keys():
			if self.gender[res] == "M" and self.chromosome == "X" and self.genotypes[res] == "AB":
				self.resultXH.append(res)
				self.result.append(res)
			if self.gender[res] == "M" and self.chromosome == "Y" and self.genotypes[res] == "AB":
				self.resultYH.append(res)
				self.result.append(res)
			if self.gender[res] == "F" and self.chromosome == "Y" and res not in self.result:
				self.resultYF.append(res)
				self.result.append(res)
			if self.chromosome == "MT" and self.genotypes[res] == "AB":
				self.resultMT.append(res)
				self.result.append(res)

		for barcode in self.result:
			if barcode in self.real_genotypes.keys():
				del self.real_genotypes[barcode]

		for barcode in self.genotypes.keys():
			if barcode in self.real_genotypes.keys():
				if self.genotypes[barcode] == 'AA' and (self.real_genotypes[barcode][1] != self.real_genotypes[barcode][0]):
					del self.real_genotypes[barcode]

				elif self.genotypes[barcode] == 'AB' and (self.real_genotypes[barcode][1] == self.real_genotypes[barcode][0]):
					del self.real_genotypes[barcode]

		if self.chromosome != 'MT' and self.chromosome != 'Y' and self.chromosome != 'X' and len(self.real_genotypes) > 0: 
			self.HW = self.Hardy_Weinberg_Equilibrium()

		if self.HW == False:
			for l in labels:
				self.result.append(l)


		self.result = list(set(self.result))



	def get_neighbours(self, mas, current, R):
		temp = []
		for i in range(len(mas)):
			angle_rad = np.arctan2(mas[i][1], mas[i][0])
			d = abs(max(degrees(angle_rad),current) - min(degrees(angle_rad),current))
			if (d <= R):
				temp.append(mas[i])
		return temp


	def get_center(self, neighbours):
		sumdegree = 0.0
		for i in range(len(neighbours)):
			sumdegree = sumdegree + degrees(np.arctan2(neighbours[i][1], neighbours[i][0]))
		center = sumdegree / len(neighbours)
		return center

	def delete_objs(self, array, neighbours):
		temp = []
		temp2 = []
		while(len(neighbours) != 0):
			temp = self.delete_obj(array, neighbours)
			array = temp[0]
			neighbours = temp[1]
		return array


	def delete_obj(self,array, neighbours):
		temp = []
		for i in range(len(array)):
			if(array[i][0] == neighbours[0][0] and array[i][1] == neighbours[0][1]):
				temp2 = array[:i] + array[i+1:]
				del neighbours[0]
				temp.append(temp2)
				temp.append(neighbours)
				break
		return temp


	def Hardy_Weinberg_Equilibrium(self):
		if len(list(self.real_genotypes.keys())) < 30:
			return True 

		alleles = { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0, 'D' : 0, 'I' : 0 }
		for key in self.real_genotypes.keys():
			alleles[self.real_genotypes[key][0]] = alleles[self.real_genotypes[key][0]] + 1
			alleles[self.real_genotypes[key][1]] = alleles[self.real_genotypes[key][1]] + 1

		if (alleles['D'] != 0 or alleles['I'] != 0) and (alleles['A'] != 0 or alleles['T'] != 0 or alleles['C'] != 0 or alleles['G'] != 0):
			return False
		if len([alleles[key] for key in alleles.keys() if alleles[key] > 0]) > 2:
			return False

		alleles_count = sum(alleles.values())
		observed_AF = {}

		for key in alleles:
			observed_AF[key] = alleles[key] / alleles_count

		expected_GF = {}
		for key in alleles.keys():
			for ke in alleles.keys():
				if ke != key:
					expected_GF[''.join(sorted(key + ke))] = 2 * observed_AF[key] * observed_AF[ke]
				elif ke == key:
					expected_GF[''.join(sorted(key + ke))] = observed_AF[key] ** 2

		observed_genotypes = {}

		for key in self.real_genotypes.keys():
			if ''.join(sorted(self.real_genotypes[key])) not in observed_genotypes.keys():
				observed_genotypes[''.join(sorted(self.real_genotypes[key]))] = 1
			else:
				observed_genotypes[''.join(sorted(self.real_genotypes[key]))] = 1 + observed_genotypes[''.join(sorted(self.real_genotypes[key]))]

		genotypes_count = sum(observed_genotypes.values())
		observed_GF = {}
		for key in observed_genotypes.keys():
			observed_GF[key] = observed_genotypes[key] / genotypes_count

		if len([observed_GF[key] for key in observed_GF.keys() if observed_GF[key] > 0]) > 3:
			return False


		sum_chi = 0
		for key in observed_GF.keys():
			sum_chi = sum_chi + ((observed_genotypes[key] - (expected_GF[key] * (alleles_count / 2)) ) ** 2 ) / (expected_GF[key] * (alleles_count / 2))

		if sum_chi > 3.85:
			return False
		else:
			return True



class QualityControl(object):
	
	def __init__(self, data, pics, snps=None, gender=None, angle=None):
		self.data = data
		self.snps = snps
		self.gender = gender
		self.angle = angle
		self.pics = pics
		self.main()


	def main(self):

		os.system("sort --field-separator='\t' --key=5 " + self.data + " > " + self.data + "_sorted")

		main_iterator = groupby(csv.reader(open(self.data + "_sorted"), delimiter='\t'),
				lambda row: row[0] if len(row) > 20 and row[27] != 'X' else None)
		with open('Full_report.txt', 'w') as report:
			with open('GT_report.txt', 'w') as gtreport:
				writer = csv.writer(report, delimiter = '\t')
				writer_gt = csv.writer(gtreport, delimiter = '\t')
				writer.writerow(['illumina_name', 'sample', 'code'])
				writer_gt.writerow(['illumina_name', 'sample', 'new GT', 'old GT'])
				for key, rows in main_iterator:
						if key != None and key != 'SNP Name' and (key in self.snps or self.snps == None):
							party_data = {}
							party_data_raw = {}
							flags = {}
							party_data_genotype = {}
							party_data_real_genotype = {}
							party_data_chromosome = {}
							notPassed = {}
							notPassedR = {}
							notPassedYF = {}
							notPassedYH = {}
							notPassedXH = {}
							HWstatus = {}
							genotypes = {}
							notPassedMT = {}
							notPassedC = {}
							notPassedIU = {}
							notPassedRaw = {}
							rows_ = [i for i in rows]
							for row in rows_:
								if len(row) > 20 and row[27] !='X' and row[27] != 'NaN' and row[28] != 'NaN':
									if key not in party_data.keys():
										party_data[row[0]] = []
										notPassed[key] = []
										notPassedIU[key] = []
										party_data_raw[row[0]] = []
										party_data_genotype[row[0]] = {}
										party_data_real_genotype[row[0]] = {}
										party_data_chromosome[row[0]] = row[16]
										if (row[2] == row[3] and row[2] != '-'):
											party_data_raw[row[0]].append({row[1]: [int(row[29]), int(row[30])] })
											party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
											party_data_genotype[row[0]][row[1]] = 'AA'
											party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
										elif row[2] != row[3] and row[2] != '-' and row[3] != '-':
											party_data_raw[row[0]].append({row[1]: [int(row[29]), int(row[30])] })
											party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
											party_data_genotype[row[0]][row[1]] = 'AB'
											party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
										elif row[2] == '-' or row[3] == '-':
											notPassed[row[0]].append(row[1])
											notPassedIU[row[0]].append(row[1])
									else:
										if row[2] == row[3] and row[2] != '-':
											party_data_raw[row[0]].append({row[1]: [int(row[29]), int(row[30])] })
											party_data_genotype[row[0]][row[1]] = 'AA'
											party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
											party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
										elif row[2] != row[3] and row[2] != '-' and row[3] != '-':
											party_data_raw[row[0]].append({row[1]: [int(row[29]), int(row[30])] })
											party_data_genotype[row[0]][row[1]] = 'AB'
											party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
											party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
										elif row[2] == '-' or row[3] == '-':
											notPassed[row[0]].append(row[1])
											notPassedIU[row[0]].append(row[1])

	
							for key in party_data.keys():
								if len(party_data[key]) > 0:
									notPassedR[key] = []
									notPassedYF[key] = []
									notPassedYH[key] = []
									notPassedXH[key] = []
									notPassedMT[key] = []
									notPassedC[key] = []
									notPassedRaw[key] = []
									genotypes[key] = {}
									obj = Forel(self.angle, party_data[key], key, party_data_genotype[key], party_data_chromosome[key], self.gender, party_data_real_genotype[key], party_data_raw[key])
	
	
									r = obj.result
									byR = obj.resultR
									byRaw = obj.resultRaw
									flags[key] =  obj.flag
									gen = obj.genotypes
									byYH = obj.resultYH
									byYF = obj.resultYF
									byXH = obj.resultXH
									byMT = obj.resultMT
									byC = obj.resultCluster
									HWstatus[key] = obj.HW
									for res in r:
										notPassed[key].append(res)
									for res in byR:
										notPassedR[key].append(res)
									for res in byYF:
										notPassedYF[key].append(res)
									for res in byYH:
										notPassedYH[key].append(res)
									for res in byXH:
										notPassedXH[key].append(res)	
									for res in byMT:
										notPassedMT[key].append(res)
									for res in byC:
										notPassedC[key].append(res)
									for res in gen.keys():
										genotypes[key][res] = gen[res]
									for res in byRaw:
										notPassedRaw[key].append(res)

									changed = []

									for row in rows_:
										if row[1] not in notPassed[row[0]] and (row[0] in genotypes.keys() and row[1] in genotypes[row[0]].keys()) and (genotypes[row[0]][row[1]] != party_data_genotype[row[0]][row[1]]):
											writer.writerow([row[0], row[1], '601'])
											writer_gt.writerow([row[0], row[1], party_data_genotype[row[0]][row[1]], genotypes[row[0]][row[1]] ])
											changed.append(row[1])
										elif row[1] not in notPassed[row[0]] and len(party_data[row[0]]) * 0.5 < len(notPassed[key]):
											writer.writerow([row[0], row[1], '602'])
										else:
											if HWstatus[row[0]] == False:
												writer.writerow([row[0], row[1], '301'])
											elif row[1] in notPassedR[row[0]] and row[1] in notPassedRaw[row[0]]:
												writer.writerow([row[0], row[1], '403'])
											elif row[1] in notPassedR[row[0]] and row[1] not in notPassedRaw[row[0]]:
												writer.writerow([row[0], row[1], '401'])
											elif row[1] in notPassedRaw[row[0]] and row[1] not in notPassedR[row[0]]:
												writer.writerow([row[0], row[1], '404'])													
											elif row[1] in notPassedYF[row[0]]:
												writer.writerow([row[0], row[1], '502'])
											elif row[1] in notPassedYH[row[0]]:
												writer.writerow([row[0], row[1], '503'])
											elif row[1] in notPassedXH[row[0]]:
												writer.writerow([row[0], row[1], '501'])
											elif row[1] in notPassedMT[row[0]]:
												writer.writerow([row[0], row[1], '503'])
											elif row[1] in notPassedC[row[0]]:
												writer.writerow([row[0], row[1], '402'])
											elif row[1] in notPassedIU[row[0]]:
												writer.writerow([row[0], row[1], '604'])

									if key in self.pics.keys():
										self.create_image(notPassed[key], party_data[key], changed, key)


	def create_image(self, failed, rows, changed, key):

		trace_green = go.Scatter(
			x = [list(z.values())[0][0] for z in rows if (list(z.keys())[0] not in failed and list(z.keys())[0] not in changed)],
			y = [list(z.values())[0][1] for z in rows if (list(z.keys())[0] not in failed and list(z.keys())[0] not in changed)],
			marker = dict(color = 'rgba(0, 255, 0, .9)', size=5),
			text = [list(z.keys())[0] for z in rows if (list(z.keys())[0] not in failed and list(z.keys())[0] not in changed)],
			mode = 'markers',
			name = 'OK'
		)

		trace_blue = go.Scatter(
			x = [list(z.values())[0][0] for z in rows if list(z.keys())[0] in failed ],
			y = [list(z.values())[0][1] for z in rows if list(z.keys())[0] in failed ],
			marker = dict(color = 'rgba(0, 0, 255, .9)', size=5),
			text = [list(z.keys())[0] for z in rows if list(z.keys())[0] in failed ],
			mode = 'markers',
			name = 'Excluded samples by QC'
		)

		trace_red = go.Scatter(
			x = [list(z.values())[0][0] for z in rows if list(z.keys())[0] in changed ],
			y = [list(z.values())[0][1] for z in rows if list(z.keys())[0] in changed ],
			marker = dict(color = 'rgba(255, 0, 0, .9)', size=5),
			text = [list(z.keys())[0] for z in rows if list(z.keys())[0] in changed ],
			mode = 'markers',
			name = 'Samples with changed genotype class'
		)


		data = [trace_green, trace_blue, trace_red]

		ymax = 0.0
		xmax = 0.0
		for k in rows:
			if list(k.values())[0][0] > xmax:
				xmax = list(k.values())[0][0]
			if list(k.values())[0][1] > xmax:
				xmax = list(k.values())[0][1]


		layout = go.Layout(
			title = key,
			yaxis = dict(
				range=[0, max(ymax, xmax)]
			),
			xaxis = dict(
				range=[0, max(ymax, xmax)]
			), 
			showlegend=True
		)

		
		figure = go.Figure(data=data, layout=layout)
		plotly.offline.plot(figure, filename=key + '.html', auto_open=False)



if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--snps", help="file with list of snps (Illumina names) for QC, all by default")
	parser.add_argument("--gender", help="file with pairs sample - gender (M or F) for extended XY QC")
	parser.add_argument("--angle", help="set angle of clusterization, float number")
	parser.add_argument("--pics", help="file with list of snps (Illumina names) for plotly graphs, none by default")
	parser.add_argument("data", help="Input data - Genome Studio report")

	args = parser.parse_args()

	snps_list = None
	pics_list = {}
	gender = None
	if args.snps:
		snps_list = {}
		for line in open(args.snps):
			row = line.split('\n')
			snps_list[row[0]] = ''

	if args.gender:
		gender = {}
		for line in open(args.gender):
			row = line.replace('\n', '').split('\t')
			gender[row[0]] = row[1]

	if args.pics:
		pics_list = {}
		for line in open(args.pics):
			row = line.split('\n')
			pics_list[row[0]] = ''

	if not args.angle:
		QualityControl(args.data, pics_list, snps_list, gender)
	else:
		QualityControl(args.data, pics_list, snps_list, gender, float(args.angle))