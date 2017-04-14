import csv
import os
import random
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from  math import atan2
from math import atan
from math import pi
from math import degrees
import numpy as np
from sklearn import linear_model
from scipy.spatial.distance import euclidean
import plotly.plotly as py
import plotly.graph_objs as go
import plotly
import argparse


class Forel(object):

	def __init__(self, angle, X, key, refgenotypes, chromosome, gender, real_genotypes):
		self.X_dict = X
		self.chromosome = chromosome
		self.key = key
		self.result = []
		self.resultR = []
		self.resultG = []
		self.resultMT = []
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
			self.R = 8.0

		X = []
		labels = []
		for dict_ in self.X_dict:
			X.append(list(dict_.items())[0][1])
			labels.append(list(dict_.items())[0][0])

		if self.gender == None:
			self.gender = {}
			for l in labels:
				self.gender[l] = 'U'

		xmax_party = 0
		for i in range(len(X)):
			if X[i][0] > 0.200 or X[i][1] > 0.200:
				mas.append(X[i])
			else:
				self.result.append(labels[i])
				self.resultR.append(labels[i])
			if X[i][0] > xmax_party:
				xmax_party = X[i][0]

		k = 0
		while(mas != []):
			current = mas[0]
			neighbours = self.get_neighbours(mas, current, self.R)
			center = self.get_center(neighbours)
			while(center[0] != current[0] and center[1] != current[1]):
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
					angle = np.rad2deg(np.arctan2(np.array(centers[i])[1], np.array(centers[i])[0]))
					angle_center = np.rad2deg(np.arctan2(np.array(list(zip(*classes[i]))[1])[0], np.array(list(zip(*classes[i]))[0])[0]))
				else:
					angle = np.rad2deg(np.arctan2(np.array(centers[i])[1], np.array(centers[i])[0]))
					angle_center = np.rad2deg(np.arctan2( ( abs(np.array(list(zip(*classes[i]))[1])[-1] - np.array(list(zip(*classes[i]))[1])[0])), abs((lr_result[-1] - lr_result[0]))))

				cental_line_angles.append(angle_center)
				central_angles.append(angle)
				lrr.append(lr)
				xy.append(1)
			else:
				lr = linear_model.LinearRegression().fit(np.array(list(zip(*classes[i]))[0]).reshape(-1, 1), np.array(list(zip(*classes[i]))[1]))
				lr_result = lr.predict(np.array(list(zip(*classes[i]))[0]).reshape(-1, 1))
				if (len(classes[i]) == 1):
					angle = np.rad2deg(np.arctan2(np.array(centers[i])[1], np.array(centers[i])[0]))
					angle_center = np.rad2deg(np.arctan2(np.array(list(zip(*classes[i]))[1])[0], np.array(list(zip(*classes[i]))[0])[0]))
				else:
					angle = np.rad2deg(np.arctan2(np.array(centers[i])[1], np.array(centers[i])[0]))
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
			elif (heterozgosity_cluster[i] == 1 and central_angles[i] > 34 and central_angles[i] < 55):
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
				self.resultG.append(res)
				self.result.append(res)
			if self.gender[res] == "M" and self.chromosome == "Y" and self.genotypes[res] == "AB":
				self.resultG.append(res)
				self.result.append(res)
			if self.gender[res] == "F" and self.chromosome == "Y" and res not in self.result:
				self.resultG.append(res)
				self.result.append(res)
			if self.chromosome == "MT" and self.genotypes[res] == "AB":
				self.resultMT.append(res)
				self.result.append(res)

		self.HW == self.Hardy_Weinberg_Equilibrium()

		if self.HW == False:
			self.flag = 3
			for l in labels:
				self.result.append(l)


		self.result = list(set(self.result))
		
		for i in range(len(X)):
			if labels[i] in self.genotypes.keys():
				if self.genotypes[labels[i]] == self.refgenotypes[labels[i]] and labels[i] not in self.result:
					cl.append(1)
				elif labels[i] in self.result:
					cl.append(2)
				else:
					cl.append(0)
			else:
				cl.append(2)


		trace_red = go.Scatter(
			x = [z[0] for z in X if cl[X.index(z)] == 0],
			y = [z[1] for z in X if cl[X.index(z)] == 0],
			text = [z + ' ' + self.gender[z] for z in labels if cl[labels.index(z)] == 0],
			marker = dict(color = 'rgba(255, 0, 0, .9)', size=5), 
			mode = 'markers',
			name = 'Samples with changed genotype class'
		)

		trace_green = go.Scatter(
			x = [z[0] for z in X if cl[X.index(z)] == 1],
			y = [z[1] for z in X if cl[X.index(z)] == 1],
			text = [z + ' ' + self.gender[z] for z in labels if cl[labels.index(z)] == 1],
			marker = dict(color = 'rgba(0, 255, 0, .9)', size=5),
			mode = 'markers',
			name = 'OK'
		)

		trace_blue = go.Scatter(
			x = [z[0] for z in X if cl[X.index(z)] == 2],
			y = [z[1] for z in X if cl[X.index(z)] == 2],
			text = [z + ' ' + self.gender[z] for z in labels if cl[labels.index(z)] == 2],
			marker = dict(color = 'rgba(0, 0, 255, .9)', size=5),
			mode = 'markers',
			name = 'Excluded samples by QC'
		)


		data = [trace_red, trace_green, trace_blue]


		shapes = []
		for i in range(len(classes)):
			if xy[i] == 0 and len(classes[i]) != 0:
				current_lr = list(lrr[i].predict(np.array(list(zip(*classes[i]))[0]).reshape(-1, 1)))
				current_cl = np.array(list(zip(*classes[i]))[0])
			   
				shapes.append( {'type' : 'line', 'x0':current_cl[0], 'x1':current_cl[-1], 
				'y0':current_lr[0], 'y1': current_lr[len(current_lr) - 1], 'line': {
				'color': 'rgb(255, 0, 0)',
				'width': 3,
			}, })
			elif xy[i] == 1 and len(classes[i]) != 0:
				current_lr = list(lrr[i].predict(np.array(list(zip(*classes[i]))[1]).reshape(-1, 1)))
				current_cl = np.array(list(zip(*classes[i]))[1])
				shapes.append( {'type' : 'line', 'x0': current_lr[0], 'x1':current_lr[-1], 
				'y0':np.array(list(zip(*classes[i]))[1])[0], 'y1':np.array(list(zip(*classes[i]))[1])[-1], 'line': {
				'color': 'rgb(255, 0, 0)',
				'width': 3,
			}, } )


		layout = go.Layout(
			title = self.key + ' chromosome ' + self.chromosome,
			yaxis = dict(
				range=[0, max(ymax, xmax)]
			),
			xaxis = dict(
				range=[0, max(ymax, xmax)]
			),
			shapes = [i for i in shapes]
		)

		
		figure = go.Figure(data=data, layout=layout)
		plotly.offline.plot(figure, filename=self.key + '.html', auto_open=False)


	def get_neighbours(self, mas, current, R):
		temp = []
		for i in range(len(mas)):
			angle_rad = np.arctan2(mas[i][1], mas[i][0])
			angle_rad2 = np.arctan2(current[1], current[0])
			d = abs(max(degrees(angle_rad),degrees(angle_rad2)) - min(degrees(angle_rad),degrees(angle_rad2)))
			if (d <= R):
				temp.append(mas[i])
		return temp


	def get_center(self, neighbours):
		sumx = 0.0
		sumy = 0.0
		for i in range(len(neighbours)):
			sumx = sumx + neighbours[i][0]
			sumy = sumy + neighbours[i][1]
		x = sumx / len(neighbours)
		y = sumy / len(neighbours)
		center = [x,y]
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
	
	def __init__(self, data, snps=None, gender=None, angle=None):
		self.data = data
		self.snps = snps
		self.gender = gender
		self.angle = angle
		self.main()

	def main(self):

		party_data = {}
		flags = {}
		party_data_genotype = {}
		party_data_real_genotype = {}
		party_data_chromosome = {}
		notPassed = {}
		notPassedR = {}
		notPassedG = {}
		HWstatus = {}
		genotypes = {}
		notPassedMT = {}
		notPassedC = {}


		with open(args.data, 'r') as f:
			reader = csv.reader(f, delimiter = '\t')
			for row in reader:
				if len(row) > 20 and row[27] !='X' and row[27] != 'NaN' and row[28] != 'NaN':
					if self.snps == None or row[0] in self.snps.keys():
						if row[0] not in party_data.keys():
							party_data[row[0]] = []
							party_data_genotype[row[0]] = {}
							party_data_real_genotype[row[0]] = {}
							party_data_chromosome[row[0]] = row[16]
							if (row[2] == row[3] and row[2] != '-'):
								party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
								party_data_genotype[row[0]][row[1]] = 'AA'
								party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
							elif row[2] != row[3] and row[2] != '-' and row[3] != '-':
								party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
								party_data_genotype[row[0]][row[1]] = 'AB'
								party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
						else:
							if row[2] == row[3] and row[2] != '-':
								party_data_genotype[row[0]][row[1]] = 'AA'
								party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
								party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })
							elif row[2] != row[3] and row[2] != '-' and row[3] != '-':
								party_data_genotype[row[0]][row[1]] = 'AB'
								party_data_real_genotype[row[0]][row[1]] = row[2] + row[3]
								party_data[row[0]].append({row[1] : [float(row[27]), float(row[28])] })


		for key in party_data.keys():
			if len(party_data[key]) > 0 and key not in notPassed.keys():
				notPassed[key] = []
				notPassedR[key] = []
				notPassedG[key] = []
				notPassedMT[key] = []
				notPassedC[key] = []
				genotypes[key] = {}
				obj = Forel(self.angle, party_data[key], key, party_data_genotype[key], party_data_chromosome[key], self.gender, party_data_real_genotype[key] )
				r = obj.result
				byR = obj.resultR
				flags[key] =  obj.flag
				gen = obj.genotypes
				byG = obj.resultG
				byMT = obj.resultMT
				byC = obj.resultCluster
				HWstatus[key] = obj.HW
				for res in r:
					notPassed[key].append(res)
				for res in byR:
					notPassedR[key].append(res)
				for res in byG:
					notPassedG[key].append(res)
				for res in byMT:
					notPassedMT[key].append(res)
				for res in byC:
					notPassedC[key].append(res)
				for res in gen.keys():
					genotypes[key][res] = gen[res]


		for key in party_data.keys():
			if len(party_data[key]) > 0:
				notPassed[key] = []
				notPassedR[key] = []
				notPassedG[key] = []
				genotypes[key] = {}
				obj = Forel(self.angle, party_data[key], key, party_data_genotype[key], party_data_chromosome[key], self.gender, party_data_real_genotype[key] )
				r = obj.result
				byR = obj.resultR
				flags[key] =  obj.flag
				gen = obj.genotypes
				byG = obj.resultG
				HWstatus[key] = obj.HW
				for res in r:
					notPassed[key].append(res)
				for res in byR:
					notPassedR[key].append(res)
				for res in byG:
					notPassedG[key].append(res)
				for res in gen.keys():
					genotypes[key][res] = gen[res]



		with open('reportSNPQUALITYSamples.txt', 'w') as report:
			for key in notPassed.keys():
				if HWstatus == False:
					for barcode in notPassed[key]:
						report.write(key + '\t' + barcode + '\t3\n')
				else:
					for barcode in notPassed[key]:
						if barcode in notPassedR[key]:
							report.write(key + '\t' + barcode + '\t0\n')
						elif barcode in notPassedG[key]:
							report.write(key + '\t' + barcode + '\t2\n')
						else:
							report.write(key + '\t' + barcode + '\t1\n')
		

		with open('reportSNPQUALITYOverall.txt', 'w') as report:     
			if self.gender == None:
				report.write('illumina_name' + '\t'+ '№ of samples' + '\t'+ 'Excluded by QC' + '\t'+ 'Excluded by R(Luminance)' + '\t'+ 'Excluded by clusterization' + '\t'+ 'flag' + '\n')
				for key in party_data.keys():
					if len(party_data[key]) > 0:
						report.write(key + '\t' + str(len(party_data[key])) + '\t' + str(len(notPassed[key])) + '\t' + str(len(notPassedR[key])) + '\t' + str(len(notPassed[key]) - len(notPassedR[key])) + '\t' + str(flags[key]) +  '\n')
			else:
				report.write('illumina_name' + '\t'+ '№ of samples' + '\t'+ 'Excluded by QC' + '\t'+ 'Excluded by R(Luminance)' + '\t' + 'Excluded by clusterization' + '\t' + 'Excluded from X\Y' + '\t' + 'flag' + '\n')
				for key in party_data.keys():
					if len(party_data[key]) > 0:
						report.write(key + '\t' + str(len(party_data[key])) + '\t' + str(len(notPassed[key])) + '\t' + str(len(notPassedR[key])) + '\t' + str(len(notPassed[key]) - len(notPassedR[key]) - len(notPassedG[key])) + '\t' + str(len(notPassedG[key])) + '\t' + str(flags[key]) +  '\n')



		with open('reportSNPQUALITYClusters.txt', 'w') as rep:
			rep.write('illumina_name' + '\t'+ 'sample' + '\t'+ 'new GT'+ '\t'+ 'old GT' + '\n')
			for key in genotypes.keys():
				for sample in genotypes[key].keys():
					if genotypes[key][sample] != party_data_genotype[key][sample]:
						rep.write(key + '\t'+ sample + '\t'+ genotypes[key][sample] + '\t'+ party_data_genotype[key][sample] + '\n')




if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--snps", help="file with list of snps (Illumina names) for QC, all by default")
	parser.add_argument("--gender", help="file with pairs sample - gender (M or F) for extended XY QC")
	parser.add_argument("--angle", help="set angle of clusterization, float number")
	parser.add_argument("data", help="Input data - Genome Studio report")
	args = parser.parse_args()

	snps_list = None
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

	if not args.angle:
		QualityControl(args.data, snps_list, gender)
	else:
		QualityControl(args.data, snps_list, gender, float(args.angle))
