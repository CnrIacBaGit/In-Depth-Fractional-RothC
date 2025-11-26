#!/usr/bin/python3
import subprocess
import numpy.random as rnd
import sys
import json
import os
# on the base of the code from https://machinelearningmastery.com/simple-genetic-algorithm-from-scratch-in-python/

# tournament selection
def selection(pop, scores, k=3):
	# first random selection
	selection_ix = rnd.randint(len(pop))
	for ix in rnd.randint(0, len(pop), k-1):
		# check if better (e.g. perform a tournament)
		if scores[ix] < scores[selection_ix] or ((scores[ix]==scores[selection_ix]) and ((pop[ix][0]+pop[ix][1])<(pop[selection_ix][0]+pop[selection_ix][1]))):
			selection_ix = ix
	return pop[selection_ix]

# crossover two parents to create two children
def crossover(p1, p2, r_cross):
	# children are copies of parents by default
	c1, c2 = p1.copy(), p2.copy()
	# check for recombination
	if rnd.rand() < float(r_cross):
		# linear combination
		r = rnd.rand()
		for i in range(len(c1)):
			if (i % len(params["mins"]))!=alt_idx:
				c1[i] = r*p1[i] + (1.0 - r)*p2[i]
				c2[i] = r*p2[i] + (1.0 - r)*p1[i]
				if params["types"][i % len(params["mins"])]==0: # integer parameter
					c1[i]=round(c1[i])
					c2[i]=round(c2[i])
	return [c1, c2]

# mutation operator
def mutation(p1, r_mut, mins, maxs):
	for i in range(len(p1)):
		if (i % len(params["mins"]))!=alt_idx:
			# check for a mutation
			if rnd.rand() < float(r_mut):
				# set to random value
				p1[i] = mins[i% len(params["mins"])]+rnd.rand()*(maxs[i% len(params["mins"])]-mins[i% len(params["mins"])])
				if params["types"][i % len(params["mins"])]==0: # integer parameter
					p1[i]=round(p1[i])

# genetic algorithm
def genetic_algorithm(objective_calc,n_bits, n_iter, n_pop, r_cross, r_mut, mins, maxs):
	# initial population
	pop = []
	scores = []
	for i in range(n_pop):
		pop.append([])
		scores.append(0)
		for j in range(n_bits):
			if i==0:
				if j==0:
					pop[i].append(float(mins[j % len(params["mins"])]))
				else:
					pop[i].append(float(maxs[j % len(params["mins"])]))
			else:
				if i==1:
					if j==0:
						pop[i].append(float(maxs[j % len(params["mins"])]))
					else:
						pop[i].append(float(mins[j % len(params["mins"])]))
				else:
					pop[i].append(float(mins[j % len(params["mins"])])+rnd.rand()*(float(maxs[j % len(params["mins"])])-float(mins[j % len(params["mins"])])))
			if params["types"][j % len(params["mins"])]==0: # integer parameter
				pop[i][j]=round(pop[i][j])
			# fix alternated values
			if (j % len(params["mins"])) == alt_idx:
				pop[i][j]=int(j/len(params["mins"]))+1 
	# calculate basic co2 and ch4 emissions
	p_min=[]
	for j in range(n_bits):
		if params["types"][j % len(params["mins"])]==0:
			p_min.append(0)
		else:
			p_min.append(float(mins[j % len(params["mins"])]))
	p_min[alt_idx]=1
	p_min[10]=1.0
	s=params["one_p"]
	params["one_p"]=0
	s0=objective_calc(p_min,0,0,0)
	params["one_p"]=s
	zero_co2=s0[1]
	zero_ch4=s0[2]
	max_co2=s0[1]-float(params["gf_target_co2"])
	max_ch4=s0[2]-float(params["gf_target_ch4"])
	# normalize goal function coefficients
	print("initial coefficients - "+str(round(float(params["gf_k1"]),3))+" "+str(round(float(params["gf_k2"]),3))+" "+str(round(float(params["gf_k3"]),3)))
	params["gf_k1"]=float(params["gf_k1"])/float(params["gf_max_cost"])
	sumk=float(params["gf_k1"])+float(params["gf_k2"])+float(params["gf_k3"])
	params["gf_k3"]=float(params["gf_k3"])/sumk
	params["gf_k2"]=float(params["gf_k2"])/sumk
	params["gf_k1"]=float(params["gf_k1"])/sumk
	if max_ch4>0:
		params["gf_k3"]=float(params["gf_k3"])/max_ch4
		params["gf_k2"]=float(params["gf_k2"])/max_ch4
		params["gf_k1"]=float(params["gf_k1"])/max_ch4
	print("max price - "+str(round(float(params["gf_max_cost"]),3))+" initial emissions above target - "+str(round(max_co2,3))+" "+str(round(max_ch4,3)))
	print("normalized coefficients - "+str(round(params["gf_k1"],5))+" "+str(round(params["gf_k2"],5))+" "+str(round(params["gf_k3"],5)))
	if int(params["one_p"])!=0: # calculate only one candidate
		ai=0
		bitmask=int(params["bitmask"])
		for j in range(n_bits):
			if params["types"][j % len(params["mins"])]!=0:
				if (bitmask & (1<<(j % len(params["mins"]))))!=0:
					pop[0][j]=float(params["one_p_vals"][ai])
					ai=ai+1
		scores[0] = objective_calc(pop[0],i,zero_co2,zero_ch4)[0]
		quit()
	# enumerate generations
	for gen in range(n_iter):
		# evaluate all candidates in the population
		for i in range(len(pop)):
			scores[i] = objective_calc(pop[i],i,zero_co2,zero_ch4)[0]
		if gen == 0:
			best = 0
			best_eval = scores[0]
		# check for new best solution
		for i in range(n_pop):
			if scores[i] < best_eval:
				best, best_eval = pop[i], scores[i]
		print("%d best f(%s) = %.3f" % (gen,  best, best_eval))
		# select parents
		selected = [selection(pop, scores) for _ in range(n_pop)]
		# create the next generation
		children = list()
		for i in range(0, n_pop, 2):
			# get selected parents in pairs
			p1, p2 = selected[i], selected[i+1]
			# crossover and mutation
			for c in crossover(p1, p2, r_cross):
				# mutation
				mutation(c, r_mut, mins, maxs)
				# store for next generation
				children.append(c)
		# replace population
		pop = children
	return [best, best_eval]

# alters input file adding (op=0) or multiplying (op=1) to specified columns col starting from first_line
def alter_input(fin,fon,op,cols,val,first_line):
	fi=open(fin,"rt")
	fo=open(fon,"wt")
	il=fi.readlines()
	# do operation
	i=0
	for l in il:
		if i>=first_line:
			l=l.split()
			for col in cols:
				l[col]=float(l[col])
				if op==0:
					l[col]=l[col]+val
				if op==1:
					l[col]=l[col]*val
				l[col]=str(l[col])
			l=" ".join(l)+"\n"
		fo.write(l)
		i=i+1
	fi.close()
	fo.close()

# price function
def calc_price(p):
	price=0
	for v in range(n_alt):
		mi=0
		for i in range(len(params["types"])):
			if params["types"][i]==0:
				if p[v*len(params["mins"])+i]!=0:
					price=price+price_fs[mi][p[v*len(params["mins"])+i]-1](p[v*len(params["mins"]):])*p[v*len(params["mins"])+10]
				mi=mi+1
	return price
# objective function
def objective_calc(pp,i,zero_co2,zero_ch4):
	print("starting "+str(pp)+" "+str(i))
	total_co2=0
	total_ch4=0
	suma=0
	err=0
	for v in range(n_alt):
		p=pp[v*len(params["mins"]):(v+1)*len(params["mins"])]
		print("starting variant "+str(p)+" "+str(v))
		# default inputs
		vgm="vgm.txt"
		et_file="ET.txt"
		prec_file="precipitation.txt"
		gw_file="gw_depth.txt"
		Ta_file="Ta.txt"
		SoC_file="SOC_inputs.txt"
		Kc_file="Kc.txt"
		lb=float(params["LB"])
		# alter inputs according to the candidate's variable values
		if int(p[0])==1: # soil removal
			alter_input(vgm,"vgm"+str(i)+".txt",0,[4,5],-float(p[1]),1)
			vgm="vgm"+str(i)+".txt"
			lb=lb-float(p[1])
		if int(p[2])==1: # groundwater rise
			alter_input(gw_file,"gw_depth"+str(i)+".txt",0,[1],float(p[3]),0)
			gw_file="gw_depth"+str(i)+".txt"
		ie=5
		if int(p[4])==1: # plantings - grass
			alter_input(Kc_file,"Kc"+str(i)+".txt",1,[1],1.0-0.4*float(p[5]),0)
			Kc_file="Kc"+str(i)+".txt"
			alter_input(SoC_file,"SOC_inputs"+str(i)+".txt",1,[1],float(p[6])*float(p[5]),0)
			SoC_file="SOC_inputs"+str(i)+".txt"
			ie=float(p[7])
		if int(p[4])==2: # plantings - grass
			alter_input(Kc_file,"Kc"+str(i)+".txt",1,[1],1.0-0.4*float(p[5]),0)
			Kc_file="Kc"+str(i)+".txt"
			alter_input(SoC_file,"SOC_inputs"+str(i)+".txt",1,[1],float(p[8])*float(p[5]),0)
			SoC_file="SOC_inputs"+str(i)+".txt"
			ie=float(p[9])
		# start simulation
		# create separate directory
		os.system("mkdir "+str(i)+" 2>/dev/null")
		# make symlinks
		os.system("ln -s $PWD/z_rothC "+str(i)+"/z_rothC"+" 2>/dev/null")
		os.system("ln -s $PWD/"+vgm+" "+str(i)+"/"+vgm+" 2>/dev/null")
		os.system("ln -s $PWD/"+et_file+" "+str(i)+"/"+et_file+" 2>/dev/null")
		os.system("ln -s $PWD/"+prec_file+" "+str(i)+"/"+prec_file+" 2>/dev/null")
		os.system("ln -s $PWD/"+gw_file+" "+str(i)+"/"+gw_file+" 2>/dev/null")
		os.system("ln -s $PWD/"+Ta_file+" "+str(i)+"/"+Ta_file+" 2>/dev/null")
		os.system("ln -s $PWD/"+SoC_file+" "+str(i)+"/"+SoC_file+" 2>/dev/null")
		os.system("ln -s $PWD/"+Kc_file+" "+str(i)+"/"+Kc_file+" 2>/dev/null")
		# run
		cwd=os.getcwd();
		cwd=cwd+"/"+str(i)
		lst=["./z_rothC","testing","3", "Zoutstep", "1", "C_q", "1", "corrector_type", "0" ,"Tm",str(params["fet"]), "Om","1", "LB",str(lb),
		"H0",str(params["H0"]),"NB",str(params["N"]), "ls_eps",str(params["ls_eps"]), "ls_max_iter",str(params["ls_max_iter"]), 
		"vgm", vgm, "et_file",et_file, "prec_file", prec_file, "gw_file",gw_file,"Ta_file",Ta_file,"SoC_file", SoC_file,"Kc_file",Kc_file,
		"inputs_exponent",str(ie),"C00",str(params["C00"]),"C01",str(params["C01"]),"C02",str(params["C02"]),"C03",str(params["C03"])];
		proc=subprocess.Popen(lst,stdout=subprocess.PIPE,cwd=cwd)
		try:
			proc.communicate(timeout=300)
		except subprocess.TimeoutExpired:
			proc.kill()
			proc.communicate()
			err=1
		if err==1:
			print("error while simulating")
			break
		# read CO2 and CH4 emissions
		fi=open(str(i)+"/out_C3.txt","rt")
		il=fi.readlines()
		fi.close()
		l=il[-1].split()
		co2=float(l[int(params["N"])+26])
		ch4=float(l[int(params["N"])+28])
		total_co2=total_co2+co2*float(p[10])
		total_ch4=total_ch4+ch4*float(p[10])
		suma=suma+float(p[10])
		print("co2:"+str(co2)+" ch4:"+str(ch4)+" area:"+str(p[10]))
		if int(params["one_p"])!=0: # stop on i-th alternate
			if v==(int(params["one_p"])-1):
				quit()
	ok=1
	if err==1:
		ok=0
	# add emissions for non-altered areas
	if suma<=float(params["maxs"][10]):
		total_co2=total_co2+(float(params["maxs"][10])-suma)*zero_co2
		total_ch4=total_ch4+(float(params["maxs"][10])-suma)*zero_ch4
	else:
		ok=0
	# calculate price
	price=calc_price(pp)
	if price>float(params["gf_max_cost"]):
		ok=0
	# calculate the level of reach of target emissions values
	rco2=total_co2-float(params["gf_target_co2"])
	if rco2<0:
		rco2=0
	rch4=total_ch4-float(params["gf_target_ch4"])
	if rch4<0:
		rch4=0
	gf=float(params["gf_k1"])*price+float(params["gf_k2"])*rco2+float(params["gf_k3"])*rch4
	# set penalty value if needed
	if ok==0:
		gf=100000
	# calculate final goal function value
	print(str(i)+":("+str(pp)+")->pr:"+str(round(price,3))+" co2:"+str(round(total_co2,3))+" rco2:"+str(round(rco2,3))+" ch4:"+str(round(total_ch4,3))+" rch4:"+str(round(rch4,3))+"->gf:"+str(round(gf,3)))
	return [gf,total_co2,total_ch4]

params_fname="z_rothC_genetic_params.json"
if len(sys.argv)>1:
	params_fname=str(sys.argv[1])
# read parameters json
print("reading params from "+params_fname)
jsf=open(params_fname)
params=json.load(jsf)
jsf.close()
# parse argv for <param_name> <param_value> pairs
if len(sys.argv)>2:
	for i in range(int((len(sys.argv)-2)/2)):
		s=str(sys.argv[2*i+3])
		if s[0].isalpha():
			s="\""+s+"\""
		params[str(sys.argv[2*i+2])]=json.loads(s)
		print("set "+str(sys.argv[2*i+2])+"="+s)
# process price functions
price_fs=[]
for i in range(len(params["prf"])):
	price_fs.append([])
	for j in range(len(params["prf"][i])):
		price_fs[i].append(eval(params["prf"][i][j]))
# process bitmask - set max to min if no fitting
bitmask=int(params["bitmask"])
for i in range(len(params["mins"])):
	if (bitmask & (1<<i))==0:
		params["maxs"][i]=params["mins"][i]
# alternative measures
alt_idx=int(params["alt_areas"])
n_alt=1
if params["types"][alt_idx]==0:
	n_alt=int(params["maxs"][alt_idx])
else:
	n_alt=1
# run optimization
best,best_eval=genetic_algorithm(objective_calc,n_alt*len(params["mins"]),int(params["gen_niter"]),int(params["gen_npop"]),params["gen_rcross"],params["gen_rmut"],params["mins"],params["maxs"])
print(best)
print(best_eval)