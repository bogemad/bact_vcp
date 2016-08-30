#!/usr/bin/env python

import gzip, os, sys, collections, re

class FASTQRead(collections.namedtuple('FASTQRead', 'name comment seq qual')):
	@property
	def fullname(self):
		if self.comment:
			return '%s %s' % (self.name, self.comment)
		else:
			return self.name

	def __repr__(self):
		if self.comment:
			return '@%s %s\n%s\n+\n%s\n' % (self.name, self.comment, self.seq, self.qual)
		else:
			return '@%s\n%s\n+\n%s\n' % (self.name, self.seq, self.qual)

	def subseq(self, start, end, comment=None):
		if self.comment:
			comment = '%s %s' % (self.comment, comment)

		return FASTQRead(self.name, comment, self.seq[start:end], self.qual[start:end])

	def clone(self, name=None, comment=None, seq=None, qual=None):
		n = name if name else self.name
		c = comment if comment else self.comment
		s = seq if seq else self.seq
		q = qual if qual else self.qual

		return FASTQRead(n, c, s, q)

	def write(self, out):
		out.write(repr(self))


def fastq_read_file(fileobj):
	name = fileobj.next().strip()[1:]

	spl = re.split(r'[ \t]', name, maxsplit=1)
	name = spl[0]
	if len(spl) > 1:
		comment = spl[1]
	else:
		comment = ''

	seq = fileobj.next().strip()
	fileobj.next()
	qual = fileobj.next().strip()

	return FASTQRead(name, comment, seq, qual)


class FASTQ(object):
	def __init__(self, fname=None, fileobj=None):
		self.fname = fname
		self._is_paired = None
		self._is_colorspace = None

		if fileobj:
			self.fileobj = fileobj
		elif fname:
			if fname == '-':
				self.fileobj = sys.stdin
			elif fname[-3:] == '.gz' or fname[-4:] == '.bgz':
				self.fileobj = gzip.open(os.path.expanduser(fname))
			else:
				self.fileobj = open(os.path.expanduser(fname))
		else:
			raise ValueError("Must pass either a fileobj or fname!")

	def tell(self):
		# always relative to uncompressed...
		return self.fileobj.tell()

	def seek(self, pos, whence=0):
		self.fileobj.seek(pos, whence)

	def fetch(self, quiet=False, callback=None):
		# if self.fname and not quiet:
			# eta = ETA(os.stat(self.fname).st_size, fileobj=self.fileobj)
		# else:
			# eta = None

		while True:
			try:
				read = fastq_read_file(self.fileobj)
				# if eta:
					# if callback:
						# eta.print_status(extra=callback())
					# else:
						# eta.print_status(extra=read.name)
				yield read

			except StopIteration:
				break

		# if eta:
			# eta.done()

	def close(self):
		if self.fileobj != sys.stdout:
			self.fileobj.close()

	def check_qualtype(self, num_to_check=10000):
		'''
		Checks a FASTQ file's quality score to see what encoding/scaling is used:
		Sanger, Solexa, or Illumina

		returns "Sanger", "Solexa", "Illumina", or "Unknown"
		'''

		pos = self.tell()

		# these are the differential values, unscaled from chr()
		sanger = (33, 74)  # default sanger is 0->40, but some newer illumina on this scale is 0->41
		solexa = (59, 104)
		illumina = (64, 104)

		sanger_count = 0
		solexa_count = 0
		illumina_count = 0
		unknown_count = 0

		checked = 0
		for read in self.fetch(quiet=True):
			if checked > num_to_check:
				break
			qmax = None
			qmin = None
			for q in [ord(x) for x in read.qual]:
				if qmin is None or q < qmin:
					qmin = q
				if qmax is None or q > qmax:
					qmax = q

			if sanger[0] <= qmin <= qmax <= sanger[1]:
				sanger_count += 1
			elif illumina[0] <= qmin <= qmax <= illumina[1]:
				illumina_count += 1
			elif solexa[0] <= qmin <= qmax <= solexa[1]:
				solexa_count += 1
			else:
				unknown_count += 1
			checked += 1

		self.seek(pos)

		if unknown_count > 0:
			return 'Unknown'  # We don't have any idea about at least one of these reads

		if solexa_count > 0:
			# If there are any reads that fall in the Solexa range,
			# this must be a Solexa scale file. This should be rare.
			return 'Solexa'

		if sanger_count > illumina_count:
			return 'Sanger'
		return 'Illumina'

	@property
	def is_colorspace(self):
		'''
		This works by scanning the first 10 reads that have sequences (aren't Ns
		or 4s). If there are any colorspace values, the entire file is called as
		colorspace.

		It's a bit overkill...
		'''

		if self._is_colorspace is not None:
			return self._is_colorspace

		pos = self.tell()
		self.seek(0)
		self._is_colorspace = None

		valid_basespace = "atcgATCG"
		valid_colorspace = "0123456"

		for read in self.fetch(quiet=True):
			if len(read.seq) < 2:
				continue

			for base in read.seq[1:]:  # skip the first base, in case there is a linker prefix
				if base in valid_colorspace:
					self._is_colorspace = True
					break
				elif base in valid_basespace:
					self._is_colorspace = False
					break
			if self._is_colorspace is not None:
				break

		self.seek(pos)
		return self._is_colorspace

	@property
	def pair_count(self):
		if self.is_paired:  # this actually does the calculation and is cached
			return self._is_paired
		return 1

	@property
	def is_paired(self):
		'''
		Determines if a FASTQ file has paired reads. This returns True is the file has
		paired reads with the same name in consecutive order.
		'''

		if self._is_paired is not None:
			return (self._is_paired > 1)

		pos = self.tell()
		self.seek(0)
		last_name = None
		count = 0

		for read in self.fetch(quiet=True):
			name = read.name.split()[0]
			if last_name:
				if name == last_name:
					count += 1
				else:
					self._is_paired = count
					self.seek(pos)

					return (self._is_paired > 1)
			else:
				last_name = name
				count = 1

		# if there are only 2 reads...
		self._is_paired = count
		self.seek(pos)
		return (self._is_paired > 1)


def unmerge(combined_fname, out_template, gz=False):
	outs = []
	if gz:
		outs.append(gzip.open('%s.1.fastq.gz' % out_template, 'w'))
	else:
		outs.append(open('%s.1.fastq' % out_template, 'w'))

	outidx = 1

	last_read = None
	fq = FASTQ(combined_fname)
	for read in fq.fetch():
		if last_read and last_read.name == read.name:
			outidx += 1
			if len(outs) < outidx:
				if gz:
					outs.append(gzip.open('%s.%s.fastq.gz' % (out_template, outidx), 'w'))
				else:
					outs.append(open('%s.%s.fastq' % (out_template, outidx), 'w'))
			read.write(outs[outidx - 1])
		else:
			outidx = 1
			read.write(outs[0])

		last_read = read

	fq.close()
	for out in outs:
		out.close()


def generator_fetch(generators):
	while True:
		try:
			yield [gen.next() for gen in generators]
		except StopIteration:
			return


def merge(fastqs, split_slashes=False, out=sys.stdout, quiet=False):
	for reads in generator_fetch([fq.fetch(quiet=quiet if i == 0 else False) for i, fq in enumerate(fastqs)]):
		cur_name = None
		for read in reads:
			name = read.name
			comment = read.comment

			if split_slashes and '/' in name:
				spl = name.split('/', 1)
				name = spl[0]
				if read.comment:
					comment = '/%s %s' % (spl[1], read.comment)
				else:
					comment = '/%s' % spl[1]

			if not cur_name:
				cur_name = name
			else:
				if name != cur_name:
					raise ValueError('Files are not paired! Expected: "%s", got "%s"!' % (cur_name, name))

			read.clone(name=name, comment=comment).write(out)

