# This program count the frequency of motifs from a batch of genomic regions.
# Genomic regions are BED format.
import sys
import os
import numpy as np
from Bio import SeqIO, Motif
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot
import random
import itertools
from time import time, strftime
import uuid
import getopt

# Functions definition
def search_motif(motiflist, seq, col, extend):
	"""search motif PWM from sequence list"""
	freq_list = []
	shift_list = []
	mf = Motif.read(open(motiflist),'jaspar-pfm')
	background = 0
	for sequence,control in itertools.izip(seq, col):
		hit = [(pos,score) for pos,score in mf.search_pwm(sequence,threshold=7.0)]
		scores = np.array([score for (pos,score) in hit])
		positions = np.array([pos for (pos,score) in hit])
		if extend != 0:
			dist = [abs(extend-base) if base >=0 else abs(-1*extend-base) for base in positions]
		freq_list.append(len(scores))
		shift_list += dist
		background += len([score for pos,score in mf.search_pwm(control,threshold=7.0)])
	return freq_list, background, len(mf), shift_list

def sliding_window(sequence,winSize,step=1):
	"""Returns a generator that will iterate through
	the defined chunks of input sequence.  Input sequence
	must be iterable."""

	# Verify the inputs
	try: it = iter(sequence)
	except TypeError:
		raise Exception("**ERROR** sequence must be iterable.")
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step must not be larger than winSize.")
	if winSize > len(sequence):
		raise Exception("**ERROR** winSize must not be larger than sequence length.")

	# Pre-compute number of chunks to emit
	numOfChunks = ((len(sequence)-winSize)/step)+1
	
	# Do the work
	for i in range(0,numOfChunks*step,step):
		yield sequence[i:i+winSize]

def generate_plot(motif_freq, motif_background, motif_shift, ext, win_size, step, tf, shift):
	"""Generate two figure for motif enrichment and spatial analysis"""
	# Generate the enrichment plot
	random.seed(1234567890)
	file_name = str(uuid.uuid4())
	pyplot.figure(1)
	for each in motif_freq:
		try:
			profile = [sum(i)/float(motif_background[each]) for i in sliding_window(motif_freq[each],win_size,step)]
		except ZeroDivisionError:
			profile = [sum(i) for i in sliding_window(motif_freq[each],win_size,step)]
		colors = [random.random() for i in range(3)]
	
		# Highlight the curve for the TF of input ChIP-seq
		if each == tf:
			width = 3
		else:
			width = 1
		pyplot.plot(profile, color=colors, linestyle='-', linewidth=width, label=each)
	pyplot.legend(prop={'size':7})
	pyplot.savefig(file_name+'.enrich.png',format='png')

	# Generate the spatial plot
	if ext != 0:
		pyplot.figure(2)
		for each in motif_freq:
			profile,bin = np.histogram(motif_shift[each], bins=range(0,ext,shift))
			colors = [random.random() for i in range(3)]
			if each == tf:
				width = 3
			else:
				width = 1
			pyplot.plot(profile, color=colors, linestyle='-', linewidth=width, label=each)
		pyplot.legend(prop={'size':7})
		pyplot.savefig(file_name+'.space.png',format='png')

def usage():
	"""showing help information"""
	print 'Usage:'
	print ' python NgramMotif.py -b <BEDfile> -m <MOTIFname> -r <reference_genome.fasta> [opts]'
	print 'Opts:'
	print ' -e <int>	:increase the BED entry by the same number base pairs in each direction (default 150)'
	print ' -w <int>	:window size for peak regions (default 1000)'
	print ' -s <int>	:step for sliding windows of peak regions (default 100)'
	print ' -f <int>	:unit of shift base pairs from motif start to peak summit (default 3)'
	print ' -n      	:transcription factor name of input ChIP-seq, the generated plot will highlight this motif (default None)'
	print ' -g      	:path of reference genome size (default hg19)'
	print ' -h      	:produce this menu'

if __name__ == "__main__":

	# set current path
	fullpath = os.path.realpath(__file__)
	filename = sys.argv[0].split('/')[-1]
	path = fullpath.split(filename)[0]

	# parameters parsing
	bedinput = None
	motifinput = None
	ext = 150 # default region: summit +/- 150bp
	win_size = 1000 # the size of peaks in each point in X axis
	step = 100 
	shift = 3 # motif default start site every 3bp shift from the peak summit
	tf_name = None
	reference = '/compbio/data/bwaIndex/hg19.fa'
	genome = path+'genome/human.hg19.genome'
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'b:m:e:w:s:f:n:r:g:h')
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-b':
			bedinput = a
		elif o == '-m':
			motifinput = a
		elif o == '-e':
			ext = int(a)
		elif o == '-w':
			win_size = int(a)
		elif o == '-s':
			step = int(a)
		elif o == '-f':
			shift = int(a)
		elif o == '-n':
			tf_name = a
		elif o == '-r':
			reference = a
		elif o == '-g':
			genome = a
		elif o == '-h':
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not bedinput or not motifinput:
		usage()
		sys.exit(2)
	
	# Timing progrom start
	start = time()
	print 'Ngram start: '+strftime("%Y-%m-%d %H:%M:%S")

	# Extract sequneces form BED file
	file_name = str(uuid.uuid4())
	os.system('slopBed -i '+bedinput+' -g '+genome+' -b '+str(ext)+' >'+file_name+'.bed')
	os.system('fastaFromBed -fi '+reference+' -bed '+file_name+'.bed -fo '+file_name+'.fa')
	os.system('shuffleBed -i '+file_name+'.bed -g '+genome+' >'+file_name+'_control.bed')
	os.system('fastaFromBed -fi '+reference+' -bed '+file_name+'_control.bed -fo '+file_name+'_control.fa')
	seqset = [record.seq for record in SeqIO.parse(file_name+'.fa','fasta')]
	colset = [record.seq for record in SeqIO.parse(file_name+'_control.fa','fasta')]
	os.remove(file_name+'.fa')
	os.remove(file_name+'.bed')
	os.remove(file_name+'_control.fa')
	os.remove(file_name+'_control.bed')

	# Search frequency of candidate motifs from each sequence
	motif_file = {}
	motif_freq = {}
	motif_background = {}
	motif_len = {}
	motif_shift = {}
	ifp = open(path+'transfact/namelist.txt')
	for line in ifp:
		item = line.rstrip().split()
		motif_file[item[0]] = item[1]
	ifp.close()
	ifp = open(motifinput)
	for line in ifp:
		item = line.rstrip()
		try:
			pfm = path+'transfact/'+motif_file[item]
		except KeyError:
			print item+' motif not in the database!'
		motif_freq[item], motif_background[item], motif_len[item], motif_shift[item] = search_motif(pfm, seqset, colset, ext)
	ifp.close()

	# Generating summary of motif searching
	print '-----------------------------------------------'
	print '# of motifs:', len(motif_len)
	print '# of sequences:', len(seqset)
	print '-----------------------------------------------'
	print 'name\tlength\toccurrence#\tbackground#'
	for each in motif_freq:
		print each+'\t'+str(motif_len[each])+'\t'+str(sum(motif_freq[each]))+'\t'+str(motif_background[each])
	print '-----------------------------------------------'

	# Generating the motif searches figures
	generate_plot(motif_freq, motif_background, motif_shift, ext, win_size, step, tf_name, shift)
	print 'Ngram finish: '+strftime("%Y-%m-%d %H:%M:%S")
	end = time()
	print 'Programing runing: '+str(end-start)+' seconds.'
