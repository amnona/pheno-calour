import requests
from csv import DictReader
import re
import codecs


def to_utf8(infile,outfile):
	BLOCKSIZE = 1048576  # or some other, desired size in bytes
	with codecs.open(infile, "r", "latin8") as sourceFile:
		with codecs.open(outfile, "w", "utf-8") as targetFile:
			while True:
				contents = sourceFile.read(BLOCKSIZE)
				if not contents:
					break
				targetFile.write(contents)


def import_data(filename, fprimer='GTGCCAGC[AC]GCCGCGGTAA', len=150):
	'''Import the phenotypedb tsv file and find the v4 region for each sequence

	Parameters
	----------
	filename: str

	Returns
	-------
	'''
	efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
	sdict = {}
	with open(filename) as dbf:
		reader = DictReader(dbf, delimiter='\t')
		for crow in reader:
			accession = crow['16S rDNA accession number']
			print(accession)
			accession = accession.strip(' .,')
			try:
				res = requests.get(efetch, params={'db':'nucleotide', 'rettype':'FASTA', 'id':accession})
			except:
				print('timeout')
				continue
			if res.status_code != 200:
				print('fetch failed for: %s' % accession)
				continue
			try:
				cseq = ''.join(res.content.decode("utf-8").split('\n')[1:])
				match=re.search(fprimer, cseq)
				cseq = cseq[match.end():]
				cseq = cseq[:len]
			except:
				print('error processing sequence: %s' % accession)
				continue
			if cseq in sdict:
				print('sequence %s already in sdict' % accession)
			else:
				sdict[cseq]={}
			for ck,cv in crow.items():
				if cv != '':
					sdict[cseq][ck] = cv
	return sdict


def export_data(sdict, outfilename):
	'''Export the phenotype database with sequence per entry

	Parameters
	----------
	sdict : dict of {sequence:phenotypes}
		from import_data
	outfilename: str
		name of the output tsv file
	'''
	print('saving to file %s' % outfilename)
	with open(outfilename,'w') as ofl:
		ofl.write('Sequence')
		columns = set()
		for cdict in sdict.values():
			for ck in cdict.keys():
				columns.add(ck)
		columns = list(columns)
		print('found %d columns' % len(columns))
		for ccol in columns:
			ofl.write('\t%s' % ccol)
		ofl.write('\n')
		numseqs=0
		for cseq, cdict in sdict.items():
			ofl.write('%s' % cseq)
			for ccol in columns:
				ofl.write('\t')
				if ccol in cdict:
					ofl.write('%s' % cdict[ccol])
			ofl.write('\n')
			numseqs += 1
		print('wrote %d seqs' % numseqs)
