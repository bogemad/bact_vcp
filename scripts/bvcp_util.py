#!/usr/bin/env python

from contextlib import contextmanager
import os, shutil, subprocess, sys

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)


def replace_dir(dir):
	if os.path.exists(dir):
		shutil.rmtree(dir)
	os.mkdir(dir)


def run_command(command, outputfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	while True:
		output = process.stdout.readline()
		if output == '' and process.poll() is not None:
			break
		if output:
			with open(outputfile,'a+') as outhandle:
				outhandle.write(output)
	if process.poll() == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), outputfile)
		sys.exit(1)


def run_po_command(command, outputfile, errorfile):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, error = stdout, stderr = process.communicate()
	if output:
		with open(outputfile,'w') as outhandle:
			outhandle.write(output)
	if error:
		with open(errorfile,'a+') as errorhandle:
			errorhandle.write(error)
	if process.returncode == 1:
		print "\n\n###########    Error    ###########\n\nThe following command failed to complete successfully:\n\n%s\n\nOutput of this error can be found in:\n\n%s\n\n" % (" ".join(command), errorfile)
		sys.exit(1)


def copy_to_results_dir(source, destination_dir):
	hold_dir = os.path.join(destination_dir,".copy")
	replace_dir(hold_dir)
	if os.path.isdir(source):
		shutil.copytree(source,os.path.join(hold_dir,os.path.basename(source)))
		os.rename(os.path.join(hold_dir,os.path.basename(source)), os.path.join(destination_dir,os.path.basename(source)))
		shutil.rmtree(hold_dir)
	elif os.path.isfile(source):
		shutil.copyfile(source,os.path.join(hold_dir,os.path.basename(source)))
		os.rename(os.path.join(hold_dir,os.path.basename(source)), os.path.join(destination_dir,os.path.basename(source)))
		shutil.rmtree(hold_dir)
	else:
		print "Error copying object: %s cannot determine if directory or file." % source
		shutil.rmtree(hold_dir)
		sys.exit(1)


def skip_line(inhandle,outhandle):
	line = inhandle.next()
	outhandle.write(line)


def add_text_to_line(inhandle,outhandle,text):
	line = inhandle.next().rstrip("\n")
	line += text
	line += "\n"
	outhandle.write(line)

