#!/usr/bin/env python
import sys
import subprocess
import argparse
import gzip
import bz2
import os
import os.path
import re
import shutil
import urllib2
from time import time, localtime, strftime
import json
#import sra_tools
import glob

"""
This script organizes a command line for either 
Unicycler, or canu as appropriate: 
    canu if only long reads (pacbio or nanopore), 
    Unicycler if illumina or iontorrent 
    and Unicycler for hybrid assemblies, eg Illumina plus PacBio
    or Spades if requested
It can auto-detect read types (illumina, iontorrent, pacbio, nanopore)
It can run trim_galore prior to assembling.
It can run Quast to generate assembly quality statistics.
TODO: properly handle different kinds of pacbio reads
TODO: verify that read type identification works in general
"""

DEFAULT_GENOME_SIZE = "5m"
Default_bytes_to_sample = 20000
MAX_SHORT_READ_LENGTH = 600
Read_id_sample = {}
Read_file_type = {}
Avg_read_length = {}
LOG = None # create a log file at start of main()
START_TIME = None
WORK_DIR = None
SAVE_DIR = None
DETAILS_DIR = None

def registerReads(reads, details, platform=None, interleaved=False, supercedes=None):
    """
    create an entry in details for these reads
    move read files to working directory to allow relative paths
    """
    LOG.write("registerReads( %s, platform=%s, interleaved=%s, supercedes=%s\n"%(reads, str(platform), str(interleaved), str(supercedes)))
    if reads in details["original_items"]:
        comment = "duplicate registration of reads %s"%reads
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None
    details['original_items'].append(reads)
    
    readStruct = {}
    readStruct["file"] = []
    #readStruct["path"] = []
    readStruct["problem"] = []
    readStruct["layout"] = 'na'
    readStruct["platform"] = 'na'
    readStruct['length_class'] = 'na'
    if ":" in reads or "%" in reads:
        if ":" in reads:
            delim = ["%", ":"][":" in reads] # ":" if it is, else "%"
        read1, read2 = reads.split(delim)
        readStruct["delim"] = delim
        if not os.path.exists(read1):
            comment = "file does not exist: %s"%read1
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        if not os.path.exists(read2):
            comment = "file does not exist: %s"%read2
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        dir1, file1 = os.path.split(read1)
        dir2, file2 = os.path.split(read2)
        if os.path.abspath(dir1) != WORK_DIR:
            LOG.write("symlinking %s to %s\n"%(read1, os.path.join(WORK_DIR,file1)))
            os.symlink(os.path.abspath(read1), os.path.join(WORK_DIR,file1))
        if os.path.abspath(dir2) != WORK_DIR:
            LOG.write("symlinking %s to %s\n"%(read2, os.path.join(WORK_DIR,file2)))
            os.symlink(os.path.abspath(read2), os.path.join(WORK_DIR,file2))
        if file1.endswith(".bz2"):
            uncompressed_file1 = file1[:-4]
            with open(os.path.join(WORK_DIR, uncompressed_file1), 'w') as OUT:
                with open(os.path.join(WORK_DIR, file1)) as IN:
                    OUT.write(bz2.decompress(IN.read()))
                    comment = "decompressing bz2 file %s to %s"%(file1, uncompressed_file1)
                    LOG.write(comment+"\n")
                    details["pre-assembly transformation"].append(comment)
                    file1 = uncompressed_file1
        if file2.endswith(".bz2"):
            uncompressed_file2 = file2[:-4]
            with open(os.path.join(WORK_DIR, uncompressed_file2), 'w') as OUT:
                with open(os.path.join(WORK_DIR, file2)) as IN:
                    OUT.write(bz2.decompress(IN.read()))
                    comment = "decompressing bz2 file %s to %s"%(file2, uncompressed_file2)
                    LOG.write(comment+"\n")
                    details["pre-assembly transformation"].append(comment)
                    file2 = uncompressed_file2
        readStruct["file"].append(file1)
        readStruct["file"].append(file2)
        #readStruct["path"].append(read1)
        #readStruct["path"].append(read2)
        # the "files" entry is an array1 of tuples, each tuple is a path, basename yielded by os.path.split
        registeredName = delim.join(sorted((file1, file2)))
    else:
        # no ':' or '%' delimiter, so a single file
        if not os.path.exists(reads):
            comment = "file does not exist: %s"%reads
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            return None
        if interleaved:
            readStruct["interleaved"] = True
        dir1, file1 = os.path.split(reads)
        if os.path.abspath(dir1) != WORK_DIR:
            LOG.write("symlinking %s to %s\n"%(reads, os.path.join(WORK_DIR,file1)))
            os.symlink(os.path.abspath(reads), os.path.join(WORK_DIR,file1))
        readStruct["file"].append(file1)
        #readStruct["path"].append(reads)
        registeredName = file1

    if supercedes:
        details['derived_reads'].append(registeredName)
        readStruct['supercedes'] = supercedes
        if supercedes in details['reads']:
            details['reads'][supercedes]['superceded_by'] = registeredName
            platform = details['reads'][supercedes]['platform']
            readStruct['platform'] = platform
            try: # swap item in platform[] with this new version
                index = details['platform'][platform].index(supercedes)
                details['platform'][platform][index] = registeredName
            except ValueError:
                comment = "Problem: superceded name %s not found in details_%s"%(registeredName, platform)
    else:
        if platform:
            readStruct['platform'] = platform
            if platform not in details['platform']:
                details["platform"][platform] = []
            details["platform"][platform].append(registeredName)
    if registeredName in details['reads']:
        comment = "registered name %s already in details[reads]"%registeredName
        LOG.write(comment+"\n")
        details["problem"].append(comment)

    details['reads'][registeredName]=readStruct
    if len(readStruct['file']) == 2:
        studyPairedReads(registeredName, details)
    else:
        studySingleReads(registeredName, details)
    return registeredName

def parseJsonParameters(args):
    """ Not fully implemented: should read assembly2 service json file """
    if not os.path.exists(args.params_json):
        raise Exception("cannot find json parameters file %s\n"%args.params_json)
    LOG.write("parseJsonParameters() time=%d\n"%time())
    with open(args.params_json) as json_file:
        data = json.load(json_file)
        #args.output_dir = data['output_file']
        args.min_contig_length = data['min_contig_length']
        if "paired_end_libs" in data:
            for pe in data["paired_end_libs"]:
                pass
        if "single_end_libs" in data:
            for pe in data["single_end_libs"]:
                if pe["platform"] == 'illumina':
                    pass
    return

def studyPairedReads(item, details):
    """
    Read both files. Verify read ID are paired. Determine avg read length. Update details['reads'].
    """
    func_start = time()
    LOG.write("studyPairedReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start - START_TIME))
    details['reads'][item]['layout'] = 'paired-end'
    details['reads'][item]['avg_len'] = 0
    details['reads'][item]['length_class'] = 'na'
    details['reads'][item]['num_reads'] = 0
    file1, file2 = item.split(":")
    if file1.endswith("gz"):
        F1 = gzip.open(os.path.join(WORK_DIR, file1))
        F2 = gzip.open(os.path.join(WORK_DIR, file2))
    elif file1.endswith("bz2"):
        F1 = bz2.BZ2File(os.path.join(WORK_DIR, file1))
        F2 = bz2.BZ2File(os.path.join(WORK_DIR, file2))
    else:
        F1 = open(os.path.join(WORK_DIR, file1))
        F2 = open(os.path.join(WORK_DIR, file2))

    line = F1.readline()
    sample_read_id = line.split(' ')[0]
    F1.seek(0)
    if sample_read_id.startswith('>'):
        F1.close()
        F2.close()
        studyFastaReads(file1, details)
        studyFastaReads(file2, details)
        return
    
    read_ids_paired = True
    seqLen1 = 0
    seqLen2 = 0
    totalReadLength = 0
    seqQualLenMatch = True
    maxReadLength = 0
    minReadLength = 1e6
    maxQualScore = chr(0)
    minQualScore = chr(255)
    readNumber = 0
    i = 0
    for line1 in F1:
        line2 = F2.readline()
        if not line2:
            line2 = ""
        if i % 4 == 0 and read_ids_paired:
            read_id_1 = line1.split(' ')[0] # get part up to first space, if any 
            read_id_2 = line2.split(' ')[0] # get part up to first space, if any 
            if not readNumber:
                sample_read_id = read_id_1
            if not read_id_1 == read_id_2:
                diff = findSingleDifference(read_id_1, read_id_2)
                if diff == None or sorted(read_id_1[diff[0]:diff[1]], read_id_2[diff[0]:diff[1]]) != ('1', '2'):
                    read_ids_paired = False
                    details['reads'][item]["problem"].append("id_mismatch at read %d: %s vs %s"%(readNumber+1, read_id_1, read_id_2))
        elif i % 4 == 1:
            seqLen1 = len(line1)-1
            seqLen2 = len(line2)-1
        elif i % 4 == 3:
            if seqQualLenMatch:
                if not (seqLen1 == len(line1)-1 and seqLen2 == len(line2)-1):
                    readId = [read_id_1, read_id_2][seqLen1 != len(line1)-1]
                    seqQualLenMatch = False
                    comment = "sequence and quality strings differ in length at read %d %s"%(readNumber, readId)
                    details['reads'][item]["problem"].append(comment)
                    LOG.write(comment+"\n")
            totalReadLength += seqLen1 + seqLen2
            maxReadLength = max(maxReadLength, seqLen1, seqLen2) 
            minReadLength = min(minReadLength, seqLen1, seqLen2)
            minQualScore = min(minQualScore + line1.rstrip() + line2.rstrip())
            maxQualScore = max(maxQualScore + line1.rstrip() + line2.rstrip())
            readNumber += 1
        i += 1

    F1.close()
    F2.close()

    avgReadLength = totalReadLength/(readNumber*2)
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 

    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, maxReadLength)
    details['reads'][item]['length_class'] = ["short", "long"][maxReadLength >= MAX_SHORT_READ_LENGTH]
    if maxReadLength >= MAX_SHORT_READ_LENGTH:
        comment = "paired reads appear to be long, expected short: %s"%item
        LOG.write(comment+"\n")
        details['reads'][item]['problem'].append(comment)

    LOG.write("duration of studyPairedReads was %d seconds\n"%(time() - func_start))
    return

def studySingleReads(item, details):
    func_start = time()
    LOG.write("studySingleReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(func_start)), func_start-START_TIME))
    details['reads'][item]['layout'] = 'single-end'
    details['reads'][item]['num_reads'] = 0
    if item.endswith("gz"):
        F = gzip.open(item)
    else:
        F = open(item)

    line = F.readline()
    sample_read_id = line.split(' ')[0] # get part up to first space, if any 
    F.seek(0)
    if sample_read_id.startswith(">"):
        F.close()
        studyFastaReads(item, details)
        return

    totalReadLength = 0
    seqQualLenMatch = True
    maxReadLength = 0
    minReadLength = 1e6
    maxQualScore = chr(0)
    minQualScore = chr(255)
    readNumber = 0
    i = 0
    interleaved = True
    prev_read_id = None
    for line in F:
        if i % 4 == 0:
            read_id = line.split(' ')[0] # get part up to first space, if any 
            if not sample_read_id:
                sample_read_id = read_id
            if interleaved and i % 8 == 0 and prev_read_id: # at every other sample ID check for matching prev, indicates interleaved
                if prev_read_id != sample_read_id:
                    diff = findSingleDifference(prev_read_id, sample_read_id)
                    if diff == None or sorted(prev_read_id[diff[0]:diff[1]], sample_read_id[diff[0]:diff[1]]) != ('1', '2'):
                        interleaved=False

        elif i % 4 == 1:
            seqLen = len(line)-1
        elif i % 4 == 3:
            qualLen = len(line)-1
            if seqQualLenMatch and (seqLen != qualLen):
                seqQualLenMatch = False
                comment = "sequence and quality strings differ in length at read %d %s"%(readNumber, read_id)
                details['reads'][item]["problem"].append(comment)
                LOG.write(comment+"\n")
            totalReadLength += seqLen
            maxReadLength = max(maxReadLength, seqLen) 
            minReadLength = min(minReadLength, seqLen)
            minQualScore = min(minQualScore + line.rstrip())
            maxQualScore = max(maxQualScore + line.rstrip())
            readNumber += 1
        i += 1
                
    if not readNumber:
        comment = "no reads found in %s\n"%item
        LOG.write(comment+"\n")
        return
    avgReadLength = totalReadLength/readNumber
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 
    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, maxReadLength)
    details['reads'][item]['length_class'] = ["short", "long"][maxReadLength >= MAX_SHORT_READ_LENGTH]
    if interleaved:
        details['reads'][item]['interleaved'] = True

    LOG.write("duration of studySingleReads was %d seconds\n"%(time() - func_start))
    return

def studyFastaReads(item, details):
    """
    assume format is fasta
    count reads, calc total length, mean, max, min
    """
    func_start = time()
    seq = ""
    seqLen = 0
    totalReadLength = 0
    maxReadLength = 0
    minReadLength = 1e6
    readNumber = 0
    F = open(item)
    for line in F:
        if line.startswith(">"):
            readNumber += 1
            if seq:
                seqLen = len(seq)
                totalReadLength += seqLen
                maxReadLength = max(maxReadLength, seqLen) 
                minReadLength = min(minReadLength, seqLen)
                seq = ""
            else:
                sample_read_id = line.split()[0]
        else:
            seq += line.rstrip()
    if seq:
        seqLen = len(seq)
        totalReadLength += seqLen
        maxReadLength = max(maxReadLength, seqLen) 
        minReadLength = min(minReadLength, seqLen)

    avgReadLength = totalReadLength/readNumber
    details['reads'][item]['avg_len'] = avgReadLength
    details['reads'][item]['max_read_len'] = maxReadLength
    details['reads'][item]['min_read_len'] = minReadLength
    details['reads'][item]['num_reads'] = readNumber
    details['reads'][item]['sample_read_id'] = sample_read_id 
    details['reads'][item]['platform'] = 'fasta'
    details['reads'][item]['inferred_platform'] = inferPlatform(sample_read_id, maxReadLength)
    details['reads'][item]['length_class'] = ["short", "long"][maxReadLength >= MAX_SHORT_READ_LENGTH]
    if item not in details['platform']['fasta']:
        details['platform']['fasta'].append(item)

    LOG.write("duration of studyFastaReads was %d seconds\n"%(time() - func_start))
    return


def inferPlatform(read_id, maxReadLength):
    """ 
    Analyze sample of text from read file and return one of:
    illumina, iontorrent, pacbio, nanopore, ...
    going by patterns listed here: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#platform-specific-fastq-files
    these patterns need to be refined and tested
    """
    if read_id.startswith(">"):
        return "fasta"
    if maxReadLength < MAX_SHORT_READ_LENGTH:
        # example illumina read id
        #@D00553R:173:HG53VBCXY:2:1101:1235:2074 1:N:0:ACAGTGAT
        parts = read_id.split(":")
        if len(parts) == 3:
            return "iontorrent"
        if len(parts) > 4:
            return "illumina"
        if re.match(r"@[A-Z]\S+:\d+:\S+:\d+:\d+:\d+:\d+ \S+:\S+:\S+:\S+$", read_id):
            return "illumina" # newer illumina
        if re.match(r"@\S+:\S+:\S+:\S+:\S+#\S+/\S+$", read_id):
            return "illumina" # older illumina
        if re.match(r"@[^:]+:[^:]+:[^:]+$", read_id):
            return "iontorrent" # 
        if re.match(r"@[SED]RR\d+\.\d+", read_id):
            return "illumina" # default short fastq type 
# NOTE: need to distinguish between PacBio CSS data types and pass to SPAdes appropriately
    if re.match(r"@\S+/\S+/\S+_\S+$", read_id): #@<MovieName> /<ZMW_number>/<subread-start>_<subread-end> :this is CCS Subread
        return "pacbio" # 
    if re.match(r"@\S+/\S+$", read_id): #@<MovieName>/<ZMW_number> 
        return "pacbio" # 
#@d5edc711-3388-4510-ace0-5d39d0d70e19 runid=999acb6b58d1c399244c42f88902c6e5eeb3cacf read=10 ch=446 start_time=2017-10-24T17:33:18Z
    if re.match(r"@[a-z0-9-]+\s+runid=\S+\s+read=\d+\s+ch=", read_id): #based on one example, need to test more 
        return "nanopore" # 
    return "pacbio" # default long fastq type

def trimGalore(details, threads=1):
    startTrimTime = time()
    LOG.write("\ntrimGalore() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    if "trim report" not in details:
        details["trim report"] = {}
    toRegister = {} # save trimmed reads to register after iterating dictionary to avoid error "dictionary changed size during iteration"
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'short' and not details['reads'][reads]['platform'] == 'fasta':
            command = ['trim_galore', '-j', str(threads), '-o', '.']
            if ':' in reads:
                splitReads = reads.split(":")
                #-j 4 -o testTrim --length 30 --paired SRR1395326_1_10pct.fastq SRR1395326_2_10pct.fastq
                command.extend(["--paired", splitReads[0], splitReads[1]])
                LOG.write("command: "+" ".join(command)+"\n")
                proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
                trimGaloreStderr = proc.stderr.read()
                return_code = proc.wait()
                LOG.write("return code = %d\n"%return_code)
                trimReads = re.findall(r"Writing validated paired-end read \d reads to (\S+)", trimGaloreStderr)
                LOG.write("regex for trimmed files returned %s\n"%str(trimReads))
                if not trimReads or len(trimReads) < 2:
                    comment = "trim_galore did not name trimmed reads output files in stderr"
                    LOG.write(comment+"\n")
                    details['reads'][reads]['problem'].append(comment)
                    continue
                comment = "trim_galore, input %s, output %s"%(reads, ":".join(trimReads))
                LOG.write(comment+"\n")
                details["pre-assembly transformation"].append(comment)
                toRegister[":".join(trimReads)] = reads

                trimReports = re.findall(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
                LOG.write("re.findall for trim reports returned %s\n"%str(trimReports))
                details["trim report"][reads]=[]
                for f in trimReports:
                    shutil.move(f, os.path.join(DETAILS_DIR, os.path.basename(f)))
                    details["trim report"][reads].append(f)
            else:
                command.append(reads)
                LOG.write("command: "+" ".join(command))
                proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
                trimGaloreStderr = proc.stderr.read()
                return_code = proc.wait()
                LOG.write("return code = %d\n"%return_code)
                trimReads = re.search(r"Writing final adapter and quality trimmed output to (\S+)", trimGaloreStderr)
                LOG.write("regex for trimmed files returned %s\n"%str(trimReads))
                if not trimReads:
                    comment = "trim_galore did not name trimmed reads output files in stderr"
                    LOG.write(comment+"\n")
                    details['reads'][reads]['problem'].append(comment)
                    continue
                trimReads = trimReads.group(1)
                comment = "trim_galore, input %s, output %s"%(reads, trimReads)
                LOG.write(comment+"\n")
                details["pre-assembly transformation"].append(comment)
                toRegister[trimReads] = reads

                trimReport = re.search(r"Writing report to '(.*report.txt)'", trimGaloreStderr)
                LOG.write("regex for trim report returned %s\n"%str(trimReport))
                if trimReport:
                    trimReport = trimReport.group(1)
                    details["trim report"][reads]=trimReport
                    shutil.move(trimReport, os.path.join(DETAILS_DIR, os.path.basename(trimReport)))
    for trimReads in toRegister:
        registerReads(trimReads, details, supercedes=toRegister[trimReads])

    LOG.write("trim_galore trimming completed, duration = %d seconds\n\n\n"%(time()-startTrimTime))

def sampleReads(filename, details=None):
    srf_time = time()
    LOG.write("sampleReads() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(srf_time)), srf_time-START_TIME))
    # figures out Read_file_type
    #return read_format and sample of read ids
    read_format = 'na'
    read_id_sample = []

    if filename.endswith("gz"):
        F = gzip.open(filename)
    elif filename.endswith("bz2"):
        F = bz2.BZ2File(filename)
    else:
        F = open(filename)
    text = F.read(Default_bytes_to_sample) #read X number of bytes for text sample
    F.close()

    LOG.write("  file text sample %s:\n%s\n\n"%(filename, text[0:50]))
    lines = text.split("\n")
    readLengths = []
    if len(lines) < 2:
        comment = "in sampleReads for %s: text sample (length %d) lacks at least 2 lines"%(filename, len(text))
        LOG.write(comment+"\n")
        details["problem"].append(comment)
    if lines[0].startswith("@"):
        read_format = 'fastq'
        for i, line in enumerate(lines):
            if i % 4 == 0:
                read_id_sample.append(line.split(' ')[0]) # get part up to first space, if any 
            elif i % 4 == 1:
                readLengths.append(len(line)-1)
    elif lines[0].startswith(">"):
        read_format = 'fasta'
        read_id_sample.append(lines[0].split()[0])
        seq = ""
        for line in lines:
            if line.startswith(">"):
                readLengths.append(len(seq))
                seq = ""
            else:
                seq += line.rstrip()
    avg_read_length = 0
    if readLengths:
        avg_read_length = sum(readLengths)/float(len(readLengths))
    LOG.write("read type %s, average read length %.1f\n"%(read_format, avg_read_length))
    return read_id_sample, avg_read_length

def findSingleDifference(s1, s2):
# if two strings differ in only a single contiguous region, return the start and end of region, else return None
    if len(s1) != len(s2):
        return None
    start = None
    end = None
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if c1 != c2:
            if end:
                return None
            if not start:
                start = i
        elif start and not end:
           end = i
    if start and not end:
        end = i+1
    return (start, end)

def categorize_anonymous_read_files(args, details):
    LOG.write("categorize_anonymous_read_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("  files=%s\n"%("\t".join(args.anonymous_reads)))

    nonSraFiles = []
    sraFiles = {}
    # first pull out any SRA files for special treatment
    for item in args.anonymous_reads:
        m = re.match("([SED]RR\d+)", item)
        if m:
            sra = m.group(1)
            if sra not in sraFiles:
                sraFiles[sra] = [] 
            sraFiles[sra].append(item)
        else:
            nonSraFiles.append(item)
    for sra in sraFiles:
        processSraFastqFiles(sorted(sraFiles[sra]), details)

    # now proceed with any non-sra files
    read_file_type = {}
    read_id_sample = {}
    singleFiles = []
    pairedFiles = []
    for item in nonSraFiles:
        if ":" in item or "%" in item:
            pairedFiles.append(item)
        else:
            read_id_sample[item], avg_read_length = sampleReads(item, details)
            read_file_type[item] = inferPlatform(read_id_sample[item][0], avg_read_length)
            comment = "interpreting %s type as %s"%(item, read_file_type[item])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)
            if read_file_type[item] is not None:
                singleFiles.append(item)

    # try to find paired files
    membersOfPairs = set()
    for i, filename1 in enumerate(singleFiles[:-1]):
        for filename2 in singleFiles[i+1:]:
            singleDiff = findSingleDifference(filename1, filename2)
# singleDiff will be not None if the strings match at all but one character (presumably '1' vs '2')
            if singleDiff and singleDiff[0] > 0 and singleDiff[1]-singleDiff[0] == 1:
                charBefore = filename1[singleDiff[0]-1]
                if charBefore.isdigit():
                    continue # changes in multi-digit numbers are not indicative of paired reads
                diffChars = (filename1[singleDiff[0]], filename2[singleDiff[0]])
                pair = None
                if diffChars[0] == '1' and diffChars[1] == '2':
                    pair = (filename1, filename2)
                elif diffChars[1] == '1' and diffChars[0] == '2':
                    pair = (filename2, filename1)
                if pair:
                    comment = "candidate paired files: %s  %s"%pair
                    LOG.write(comment+"\n")
                    details["problem"].append(comment)
                    if read_file_type[filename1] != read_file_type[filename2]:
                        comment = "Discordant fileTypes for %s(%s) vs %s(%s)"%(filename1, read_file_type[filename1], filename2, read_file_type[filename2])
                        LOG.write(comment+"\n")
                        details["problem"].append(comment)
                        continue
                    pairedFiles.append(pair[0] + ":" + pair[1])

    # now go over all pairs to test for matching types and matching read IDs
    valid_pairs = set()
    valid_singles = set()
    for item in pairedFiles:
        if ":" in item:
            filename1, filename2 = item.split(":") 
        elif "%" in item:
            filename1, filename2 = item.split("%") 
        else:
            comment = "failed to find ':' or '%%' in file pair: %s" % (item)
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            continue
        if filename1 not in read_file_type:
            read_id_sample[filename1], avg_read_length = sampleReads(filename1, details)
            read_file_type[filename1] = inferPlatform(read_id_sample[filename1][0], avg_read_length)
            comment = "interpreting %s type as %s"%(filename1, read_file_type[filename1])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)
        if filename2 not in read_file_type:
            read_id_sample[filename2], avg_read_length = sampleReads(filename2, details)
            read_file_type[filename2] = inferPlatform(read_id_sample[filename2][0], avg_read_length)
            comment = "interpreting %s type as %s"%(filename2, read_file_type[filename2])
            LOG.write(comment+"\n")
            details["pre-assembly transformation"].append(comment)

        read_types_match = True
        # test if read types are the same
        if read_file_type[filename1] != read_file_type[filename2]:
            comment = "Discordant fileTypes for %s(%s) vs %s(%s)"%(filename1, read_file_type[filename1], filename2, read_file_type[filename2])
            LOG.write(comment+"\n")
            details["problem"].append(comment)
            read_types_match = False
        read_file_type[item] = read_file_type[filename1] #easier to retrieve later

        # test if read IDs match between files
        ids_paired = True
        for idpair in zip(read_id_sample[filename1], read_id_sample[filename2]):
            if idpair[0] == idpair[1]:
                continue
            diff = findSingleDifference(idpair[0], idpair[1])
            # diff reports start and end of contiguous region of different characters (if only one)
            if not diff or sorted(idpair[0][diff[0]:diff[1]], idpair[1][diff[0]:diff[1]]) != ('1', '2'):
                ids_paired = False
                comment = "Read IDs do not match for %s(%s) vs %s(%s)"%(filename1, idpair[0], filename2, idpair[1])
                LOG.write(comment+"\n")
                details["problem"].append(comment)
                singleFiles.extend((filename1, filename2)) #move over to single files
                break
        if read_types_match and ids_paired:
            valid_pairs.add(item)
            membersOfPairs.add(filename1)
            membersOfPairs.add(filename2)
        else: #move over to single files
            valid_singles.add(filename1)
            valid_singles.add(filename2)
    
    # some items on singleFiles may not be valid (may be paired up)
    for item in singleFiles:
        if item not in membersOfPairs:
            valid_singles.add(item)

    for item in valid_pairs.union(valid_singles):
        registerReads(item, details, platform=read_file_type[item], interleaved = args.interleaved and item in args.interleaved)

    return

def get_sra_runinfo(run_accession, log=None):
    """ take sra run accession (like SRR123456)
    Use edirect tools esearch and efetch to get metadata (sequencing platform, etc).
    return dictionary with keys like: spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path.....
    Altered from versionin sra_tools to handle case of multiple sra runs returned by query.
    If efetch doesn't work, try scraping the web page.
    """
    LOG.write("get_sra_runinfo(%s)\n"%run_accession)
    if run_accession.endswith(".sra"):
        run_accession = run_accession[:-4]
    runinfo = None
    runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="+run_accession
    text = urllib2.urlopen(runinfo_url).read()
    if text.startswith("Run"):
        lines = text.split("\n")
        keys   = lines[0].split(",")
        for line in lines[1:]:  # there might be multiple rows, only one of which is for our sra run accession
            if line.startswith(run_accession):
                values = line.split(",")
                runinfo = dict(zip(keys, values))
    if runinfo:
        if runinfo['Platform'] not in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
            # rescue case of mis-alignment between keys and values
            if log:
                log.write("problem in get_sra_runinfo: sra.cgi returned:\n"+text+"\n")
            for val in values:
                if val in ('ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT'):
                    runinfo['Platform'] = val
                    break
        if not runinfo['LibraryLayout'].startswith(('PAIRED', 'SINGLE')):        
            if log:
                log.write("Need to search for LibraryLayout: bad value: %s\n"%runinfo['LibraryLayout'])
            for val in values:
                if val.startswith(('PAIRED', 'SINGLE')):
                    runinfo['LibraryLayout'] = val
                    break

    if not runinfo:
        if log:
            log.write("Problem, normal runinfo request failed. Trying alternative from web page.\n")
        # screen-scrape
        runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run="+run_accession
        text = urllib2.urlopen(runinfo_url).read()
        runinfo = {}
        if re.search("<td>Illumina", text, re.IGNORECASE):
            runinfo['Platform'] = 'ILLUMINA'
        elif re.search("<td>PacBio", text, re.IGNORECASE):
            runinfo['Platform'] = 'PACBIO_SMRT'
        elif re.search("<td>Oxford", text, re.IGNORECASE):
            runinfo['Platform'] = 'OXFORD_NANOPORE'
        elif re.search("<td>Ion Torrent", text, re.IGNORECASE):
            runinfo['Platform'] = 'ION_TORRENT'

        if re.search("<td>SINGLE</td>", text, re.IGNORECASE):
            runinfo['LibraryLayout'] = 'SINGLE'
        elif re.search("<td>PAIRED", text, re.IGNORECASE):
            runinfo['LibraryLayout'] = 'PAIRED'
    return runinfo

def fetch_one_sra(sra, run_info=None, log=sys.stderr, usePrefetch=False):
    """ requires run_info to know which program to use
    """
    if not run_info:
        run_info = get_sra_runinfo(sra, log)

    if usePrefetch:
        command = ["prefetch", sra]
        log.write("command = "+" ".join(command)+"\n")
        return_code = subprocess.call(command, shell=False, stderr=log)
        log.write("return_code = %d\n"%(return_code))

    command = ["fasterq-dump", "--split-files", sra] # but not appropriate for pacbio or nanopore
    if run_info['Platform'].startswith("PACBIO") or run_info['Platform'].startswith("OXFORD_NANOPORE"):
        command = ['fastq-dump', sra]
    stime = time()
    log.write("command = "+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=log)
    log.write("return_code = %d\n"%(return_code))
    if return_code != 0:
        log.write("Problem, return code was %d\n"%(return_code))

        log.write("Try one more time.\n")
        return_code = subprocess.call(command, shell=False, stderr=LOG)
        log.write("Return code on second try was %d\n"%return_code)
        if return_code != 0:
            LOG.write("Giving up on %s\n"%sra)
    LOG.write("fetch_one_sra time=%d seconds\n"%(time()-stime))
    return

def fetch_sra_files(sra_ids, details):
    """ 
    fetch each sra item and register the reads
    """
    LOG.write("fetch_sra_files() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("sra_ids="+" ".join(sra_ids)+"\n")
    for sra in sra_ids:
        sraFull = sra
        sra = sraFull.replace(".sra", "")
        runinfo = get_sra_runinfo(sra)
        if not runinfo:
            LOG.write("runinfo for %s was empty, giving up\n"%sra)
            continue
        LOG.write("Runinfo for %s reports platform = %s and LibraryLayout = %s\n"%(sra, runinfo["Platform"], runinfo['LibraryLayout']))

        fetch_one_sra(sraFull, run_info=runinfo, log=LOG)
        fastqFiles = glob.glob(sra+"*fastq")
        LOG.write("Fastq files from sra: %s\n"%str(fastqFiles))
        processSraFastqFiles(fastqFiles, details, runinfo)
    return

def processSraFastqFiles(fastqFiles, details, run_info=None):
    """ manipulate multiple (or single) fastq files from one SRA runId and register reads """
    comment = "processSraFastqFiles(%s)"%",".join(fastqFiles)
    details['problem'].append(comment)
    LOG.write(comment+"\n")
    item = None
    m = re.match(r"([SED]RR\d+)", os.path.basename(fastqFiles[0]))
    if not m:
        comment = "supposed sra fastq file does not start with [SED]RRnnnn"
        details['problem'].append(comment)
        LOG.write(comment+"\n")
        return
    sra = m.group(1)
    for fq in fastqFiles:
        if not os.path.basename(fq).startswith(sra):
            comment = "Problem: not all fastqFiles passed to processSraFastqFiles() begin with %s: %s"%(sra, ",".join(fastqFiles))
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            return
    if not run_info:
        run_info = get_sra_runinfo(sra)

    if run_info['LibraryLayout'].startswith("PAIRED"):
        if len(fastqFiles) == 2:
            item = ":".join(sorted(fastqFiles)[:2])
            comment = "runinfo[LibraryLayout] == PAIRED: item = %s"%item
            details['problem'].append(comment)
            LOG.write(comment+"\n")
        else:
            comment = "for PAIRED library %s, number of files was %s, expected 2: %s"%(sra, len(fastqFiles), str(fastqFiles))
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            if len(fastqFiles) == 1:
                item = fastqFiles[0] # interpret as single-end, perhaps an SRA metadata mistake
                comment = "interpret library %s as single-end"%sra
                details['problem'].append(comment)
                LOG.write(comment+"\n")
    if not item: # library layout single or failed above
        if len(fastqFiles) == 1:
            item = fastqFiles[0]
            comment = "runinfo[LibraryLayout] == %s: item = %s"%(run_info['LibraryLayout'], item)
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            
        elif len(fastqFiles) > 1:
            comment = "LibraryLayout=%s; Platform=%s: multiple files = %s"%(run_info['LibraryLayout'], run_info['Platform'], ",".join(fastqFiles))
            details['problem'].append(comment)
            LOG.write(comment+"\n")

            concatenate_command = "cat %s > %s/%s.fastq"%(" ".join(fastqFiles), WORK_DIR, sra)
            LOG.write("concatenate command:"+concatenate_command+"\n")
            subprocess.call(concatenate_command, shell=True)
            item = sra+".fastq"
            comment = "for library %s, list of files was %s, concatenated to %s"%(sra, str(fastqFiles), item)
            details['problem'].append(comment)
            LOG.write(comment+"\n")
        if not item:
            comment = "for %s no fastq file found"%sra
            details['problem'].append(comment)
            LOG.write(comment+"\n")
            # failed on that sra
   
    platform = None
    if run_info["Platform"] == "ILLUMINA":
        platform = "illumina"
    elif run_info["Platform"] == "ION_TORRENT":
        platform = "iontorrent"
    elif run_info["Platform"] == "PACBIO_SMRT":
        platform = "pacbio"
    elif run_info["Platform"] == "OXFORD_NANOPORE":
        platform = "nanopore"
    if not platform:
        ids, avg_length = sampleReads(fastqFiles[0], details)
        platform = ["illumina", "pacbio"][avg_length > MAX_SHORT_READ_LENGTH]
    registerReads(item, details, platform=platform)
    return

def writeSpadesYamlFile(details):
    LOG.write("writeSpadesYamlFile: elapsed seconds = %f\n"%(time()-START_TIME))
    outfileName = "spades_yaml_file.txt"
    OUT = open(outfileName, "w")
    OUT.write("[\n")
    
    for platform in details['platform']:
        LOG.write(platform+": "+", ".join(details['platform'][platform])+"\n")
    
    single_end_reads = []
    paired_end_reads = [[], []]
    mate_pair_reads = [[], []]
    interleaved_reads = []

    shortReadItems = []
    if 'illumina' in details['platform']:
        shortReadItems.extend(details['platform']['illumina'])
    if 'iontorrent' in details['platform']:
        shortReadItems.extend(details['platform']['iontorrent'])
    for item in shortReadItems:
        if ":" in item:
            f = details['reads'][item]['file'][0]
            paired_end_reads[0].append(f)
            f = details['reads'][item]['file'][1]
            paired_end_reads[1].append(f)
        elif "%" in item:
            f = details['reads'][item]['file'][0]
            mate_pair_reads[0].append(f)
            f = details['reads'][item]['file'][1]
            mate_pair_reads[1].append(f)
        else:
            f = details['reads'][item]['file'][0]
            if 'interleaved' in details['reads'][item]:
                interleaved_reads.append(f)
            else:
                single_end_reads.append(f)

    precedingElement=False
    if single_end_reads:
        OUT.write("  {\n    type: \"single\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(single_end_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if interleaved_reads:
        OUT.write("  {\n    type: \"paired-end\",\n    interlaced reads: [\n        \"")
        OUT.write("\",\n        \"".join(interleaved_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if paired_end_reads[0]:
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"fr\",\n")
        OUT.write("    type: \"paired-end\",\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(paired_end_reads[0]))
        OUT.write("\"\n    ],\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(paired_end_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement = True
    if mate_pair_reads[0]:
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    orientation: \"rf\",\n")
        OUT.write("    type: \"mate-pairs\",\n")
        OUT.write("    left reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[0]))
        OUT.write("\"\n    ]\n")
        OUT.write("    right reads: [\n        \""+"\",\n        \"".join(mate_pair_reads[1]))
        OUT.write("\"\n    ]\n")
        OUT.write("  }\n")
        precedingElement = True
    if details['platform']['pacbio']:
        pacbio_reads = []
        for item in details['platform']['pacbio']:
            f = details['reads'][item]['file'][0]
            pacbio_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"pacbio\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(pacbio_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if details['platform']['nanopore']:
        nanopore_reads = []
        for item in details['platform']['nanopore']:
            f = details['reads'][item]['file'][0]
            nanopore_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"nanopore\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(nanopore_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True
    if details['platform']['fasta']:
        fasta_reads = []
        for item in details['platform']['fasta']:
            f = details['reads'][item]['file'][0]
            fasta_reads.append(f)
        if precedingElement:
            OUT.write(",\n")
        OUT.write("  {\n    type: \"untrusted-contigs\",\n    single reads: [\n        \"")
        OUT.write("\",\n        \"".join(fasta_reads))
        OUT.write("\"\n    ]\n  }\n")
        precedingElement = True

    OUT.write("]\n")
    OUT.close()
    return(outfileName)    

def runQuast(contigsFile, args, details):
    LOG.write("runQuast() time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    quastDir = "quast_out"
    quastCommand = ["quast.py",
                    "-o", quastDir,
                    "-t", str(args.threads),
                    "--min-contig", str(args.min_contig_length),
                    contigsFile]
    LOG.write("running quast: "+" ".join(quastCommand)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(quastCommand, shell=False, stdout=FNULL)
    LOG.write("return code = %d\n"%return_code)
    if return_code == 0:
        shutil.move(os.path.join(quastDir, "report.html"), os.path.join(DETAILS_DIR, args.prefix+"quast_report.html"))
        shutil.move(os.path.join(quastDir, "report.tsv"), os.path.join(DETAILS_DIR, args.prefix+"quast_report.tsv"))
        shutil.move(os.path.join(quastDir, "report.txt"), os.path.join(DETAILS_DIR, args.prefix+"quast_report.txt"))
        shutil.move(os.path.join(quastDir, "transposed_report.txt"), os.path.join(DETAILS_DIR, args.prefix+"quast_transposed_report.txt"))
        shutil.move(os.path.join(quastDir, "transposed_report.tsv"), os.path.join(DETAILS_DIR, args.prefix+"quast_transposed_report.tsv"))
        details["quast_transposed_txt"] = "details/"+args.prefix+"quast_transposed_report.txt"
        details["quast_transposed_tsv"] = "details/"+args.prefix+"quast_transposed_report.tsv"
        details["quast_txt"] = "details/"+args.prefix+"quast_report.txt"
        details["quast_tsv"] = "details/"+args.prefix+"quast_report.tsv"
        details["quast_html"] = "details/"+args.prefix+"quast_report.html"

def filterContigsByMinLength(inputContigs, details, min_contig_length=300, min_contig_coverage=5, threads=1, prefix=""):
    """ 
    Write only sequences at or above min_length to output file.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("filterContigsByMinLength(%s) \n"%(inputContigs))
    shortReadDepth = None
    longReadDepth = None
    bamFiles = []
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'short':
            bam = runBowtie(inputContigs, reads, details, threads=threads, outformat='bam')
            if bam:
                bamFiles.append(bam)
    if bamFiles:
        shortReadDepth = calcReadDepth(bamFiles)
    bamFiles = []
    for reads in details['reads']:
        if details['reads'][reads]['length_class'] == 'long':
            bam = runMinimap(inputContigs, reads, details, threads=threads, outformat='bam')
            if bam:
                bamFiles.append(bam)
    if bamFiles:
        longReadDepth = calcReadDepth(bamFiles)
    report = []
    report.append("min_contig_length = %d"%min_contig_length)
    report.append("min_contig_coverage = %.1f"%min_contig_coverage)
    bad_contigs = []
    details['circular_contigs'] = []
    no_coverage_available = False
    if not shortReadDepth and not longReadDepth:
        no_coverage_available = True
        comment = "No read coverage information available"
        LOG.write(comment+"\n")
        report.append(comment)
    else:
        if shortReadDepth:
            comment = "Short read coverage information available"
            LOG.write(comment+"\n")
            report.append(comment)
        if longReadDepth:
            comment = "Long read coverage information available"
            LOG.write(comment+"\n")
            report.append(comment)
    suboptimalContigs = ""
    num_good_contigs = num_bad_contigs = 0
    outputContigs = re.sub(r"\..*", "_depth_cov_filtered.fasta", inputContigs)
    LOG.write("writing filtered contigs to %s\n"%outputContigs)
    with open(inputContigs) as IN:
        with open(outputContigs, 'w') as OUT:
            seqId=None
            seq = ""
            contigIndex = 1
            line = "1"
            while line:
                line = IN.readline()
                m = re.match(r">(\S+)", line)
                if m or not line: 
                    if seq:
                        contigId = ">"+prefix+"contig_%d"%contigIndex
                        contigInfo = " length %5d"%len(seq)
                        contigIndex += 1
                        contigCoverage = 0
                        if shortReadDepth and seqId in shortReadDepth:
                            meanDepth, normalizedDepth = shortReadDepth[seqId]
                            contigInfo += " coverage %.01f normalized_cov %.2f"%(meanDepth, normalizedDepth)
                            contigCoverage = meanDepth
                        if longReadDepth and seqId in longReadDepth:
                            meanDepth, normalizedDepth = longReadDepth[seqId]
                            contigInfo += " longread_coverage %.01f normalized_longread_cov %.2f"%(meanDepth, normalizedDepth)
                            contigCoverage = max(meanDepth, contigCoverage)
                        if "contigCircular" in details and contigIndex in details["contigCircular"]:
                            contigInfo += " circular=true"
                            details['circular_contigs'].append(contigId+contigInfo)
                        if len(seq) >= min_contig_length and (no_coverage_available or contigCoverage >= min_contig_coverage):
                            OUT.write(contigId+contigInfo+"\n")
                            for i in range(0, len(seq), 60):
                                OUT.write(seq[i:i+60]+"\n")
                            num_good_contigs += 1
                        else:
                            suboptimalContigs += contigId+contigInfo+"\n"+seq+"\n"
                            bad_contigs.append(contigId+contigInfo)
                            num_bad_contigs += 1
                        seq = ""
                    if m:
                        seqId = m.group(1)
                elif line:
                    seq += line.rstrip()
    report.append("Number of contigs above thresholds: %d"%num_good_contigs)
    report.append("Number of contigs below thresholds: %d"%num_bad_contigs)
    if suboptimalContigs:
        details['bad_contigs'] = bad_contigs
        suboptimalContigsFile = "contigs_below_length_coverage_threshold.fasta"
        report.append("bad contigs written to "+suboptimalContigsFile)
        details['suboptimal_contigs'] = suboptimalContigsFile
        suboptimalContigsFile = os.path.join(DETAILS_DIR, suboptimalContigsFile)
        with open(suboptimalContigsFile, "w") as SUBOPT:
            SUBOPT.write(suboptimalContigs)
    if os.path.getsize(outputContigs) < 10:
        LOG.write("failed to generate outputContigs, return None\n")
        return None
    details['contig_filtering'] = report
    comment = "trimContigsByMinLength, input %s, output %s"%(inputContigs, outputContigs)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return outputContigs

def runBandage(gfaFile, details):
    imageFormat = ".svg"
    retval = None
    if os.path.exists(gfaFile):
        plotFile = gfaFile.replace(".gfa", ".plot"+imageFormat)
        command = ["Bandage", "image", gfaFile, plotFile]
        LOG.write("Bandage command =\n"+" ".join(command)+"\n")
        try:
            return_code = subprocess.call(command, shell=False, stderr=LOG)
            LOG.write("return code = %d\n"%return_code)
            if return_code == 0:
                retval = plotFile
            else:
                LOG.write("Error creating Bandage plot\n")
        except OSError as ose:
            comment = "Problem running Bandage: "+str(ose)
            LOG.write(comment+"\n")
            details['problem'].append(comment)
    return retval

def runUnicycler(details, threads=1, min_contig_length=0, prefix=""):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("runUnicycler: elapsed seconds = %f\n"%(time()-START_TIME))
    command = ["unicycler", "-t", str(threads), "-o", '.']
    if min_contig_length:
        command.extend(("--min_fasta_length", str(min_contig_length)))
    command.extend(("--keep", "0")) # keep only assembly.gfa, assembly.fasta and unicycler.log
    command.append("--no_pilon")

    # put all read files on command line, let Unicycler figure out which type each is
    # apparently unicycler can only accept one read set in each class (I tried multiple ways to submit 2 paired-end sets, failed)
    short1 = None
    short2 = None
    unpaired = None
    long_reads = None
    for item in details['reads']:
        if 'superceded_by' in details['reads'][item]:
            continue
        files = details['reads'][item]['file'] 
        if details['reads'][item]['length_class'] == 'short':
            if len(files) > 1:
                short1 = files[0]
                short2 = files[1]
            else:
                unpaired = files[0]
        else:
            long_reads = files[0]

    if short1:
        command.extend(("--short1", short1, "--short2", short2))
    if unpaired:
        command.extend(("--unpaired", unpaired))
    if long_reads:
        command.extend(("--long", long_reads))
    # it is not quite right to send iontorrent data to spades through unicycler because the --iontorrent flag to spades will not be set

    LOG.write("Unicycler command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    unicyclerStartTime = time()
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big and unicycle.log is better
        return_code = subprocess.call(command, shell=False, stdout=FNULL)
    LOG.write("return code = %d\n"%return_code)

    unicyclerEndTime = time()
    elapsedTime = unicyclerEndTime - unicyclerStartTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)

    details["assembly"] = { 
                'assembly_elapsed_time' : elapsedHumanReadable,
                'assembler': 'unicycler',
                'command_line': " ".join(command)
                }

    LOG.write("Duration of Unicycler run was %s\n"%(elapsedHumanReadable))

    if os.path.exists("unicycler.log"):
        unicyclerLogFile = prefix+"unicycler.log"
        shutil.move("unicycler.log", os.path.join(DETAILS_DIR, unicyclerLogFile))

    if not os.path.exists("assembly.fasta"):
        comment = "Unicycler failed to generate assembly file. Check "+unicyclerLogFile
        LOG.write(comment+"\n")
        details["assembly"]["outcome"] = comment
        details["problem"].append(comment)
        return None
    details["contigCircular"] = []
    with open("assembly.fasta") as F:
        contigIndex = 1
        for line in F:
            if line.startswith(">"):
                if "circular=true" in line:
                    details["contigCircular"].append(contigIndex) 
                contigIndex += 1

    assemblyGraphFile = prefix+"assembly_graph.gfa"
    shutil.move("assembly.gfa", os.path.join(DETAILS_DIR, assemblyGraphFile))

    contigsFile = "contigs.fasta"
    shutil.move("assembly.fasta", contigsFile) #rename to canonical name
    details["assembly"]["contigs.fasta size:"] = os.path.getsize(contigsFile)
    return contigsFile

def runSpades(details, args):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write("runSpades: elapsed seconds = %f\n"%(time()-START_TIME))

    if args.illumina and args.iontorrent:
        comment = "SPAdes is not meant to process both Illumina and IonTorrent reads in the same run"
        details["problem"].append(comment)
        LOG.write(comment+"\n")
    command = ["spades.py", "--threads", str(args.threads), "-o", "."]
    if args.recipe == 'single-cell':
        command.append("--sc")
    if args.iontorrent:
        command.append("--iontorrent") # tell SPAdes that this is the read type
    yamlFile = writeSpadesYamlFile(details)
    command.extend(["--dataset", yamlFile])
    if args.trusted_contigs:
        command.extend(["--trusted-contigs", args.trusted_contigs])
    if args.untrusted_contigs:
        command.extend(["--untrusted-contigs", args.untrusted_contigs])
    if args.memory:
        command.extend(["-m", str(args.memory)])
    if args.recipe == "meta-spades":
        command.append("--meta")
        #
        # Validate arguments for metagenomic spades. It can only run with
        # a single paired-end library.
        #
        #if len(single_end_reads) > 0 or len(paired_end_reads[0]) > 1:
        #    sys.stderr.write("SPAdes in metagenomics mode can only process a single paired-end read file\n")
        #    sys.exit(1);
    if args.recipe == "plasmid-spades":
        command.append("--plasmid")
    LOG.write("SPAdes command =\n"+" ".join(command)+"\n")
    LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    spadesStartTime = time()

    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("return code = %d\n"%return_code)

    contigsFile = "contigs.fasta"
    if return_code and not os.path.exists(contigsFile):
        comment = "spades return code = %d, see if we can restart"%return_code
        LOG.write(comment+"\n")
        details['problem'].append(comment)
        #construct list of kmer-lengths to try assembling at, omitting the highest one (that may have caused failure)
        kdirs = glob.glob("K*")
        if len(kdirs) > 1:
            knums=[]
            for kdir in kdirs:
                m = re.match("K(\d+)$", kdir)
                if m:
                    k = int(m.group(1))
                    knums.append(k)
            knums = sorted(knums) 
            next_to_last_k = knums[-1]
            kstr = str(knums[0])
            for k in knums[1:-1]:
                kstr += ","+str(k)
            command.extend("--restart-from", "k%d"%next_to_last_k, "-k", kstr) 
            LOG.write("restart: SPAdes command =\n"+" ".join(command)+"\n")

            with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
                return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
            LOG.write("return code = %d\n"%return_code)

    spadesEndTime = time()
    elapsedTime = spadesEndTime - spadesStartTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)

    details["assembly"] = { 
        'assembly_elapsed_time' : elapsedHumanReadable,
        'assembler': 'spades',
        'command_line': " ".join(command)
        }

    LOG.write("Duration of SPAdes run was %s\n"%(elapsedHumanReadable))
    spadesLogFile = args.prefix+"spades.log"
    try:
        shutil.move("spades.log", os.path.join(DETAILS_DIR, spadesLogFile))
        assemblyGraphFile = args.prefix+"assembly_graph.gfa"
        shutil.move("assembly_graph_with_scaffolds.gfa", os.path.join(DETAILS_DIR, assemblyGraphFile))
    except Exception as e:
        LOG.write(str(e))
    if not os.path.exists(contigsFile):
        comment = "SPAdes failed to generate contigs file. Check "+spadesLogFile
        LOG.write(comment+"\n")
        details["assembly"]["outcome"] = comment
        details["problem"].append(comment)
        return None
    details["assembly"]["contigs.fasta size:"] = os.path.getsize(contigsFile)
    return contigsFile

def runMinimap(contigFile, longReadFastq, details, threads=1, outformat='sam'):
    #LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    """
    Map long reads to contigs by minimap2 (read paf-file, readsFile; generate paf file).
    """
    LOG.write("runMinimap(%s, %s, details, %d, %s)\n"%(contigFile, longReadFastq, threads, outformat))
    # index contig sequences
    contigIndex = contigFile.replace(".fasta", ".mmi")
    command = ["minimap2", "-t", str(threads), "-d", contigIndex, contigFile] 
    tempTime = time() 
    LOG.write("minimap2 index command:\n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("minimap2 index return code = %d, time = %d seconds\n"%(return_code, time() - tempTime))
    if return_code != 0:
        return None

    # map long reads to contigs
    contigSam = contigFile.replace(".fasta", ".sam")
    command = ["minimap2", "-t", str(threads), "-a", "-o", contigSam, contigFile, longReadFastq]
    tempTime = time()
    LOG.write("minimap2 map command:\n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stderr=FNULL)
    LOG.write("minimap2 map return code = %d, time = %d seconds\n"%(return_code, time() - tempTime))
    if return_code != 0:
        return None

    if outformat == 'sam':
        LOG.write('runMinimap returning %s\n'%contigSam)
        return contigSam

    else:
        contigsBam = convertSamToBam(contigSam, details, threads=threads)
        LOG.write('runMinimap returning %s\n'%contigsBam)
        return contigsBam
            

def convertSamToBam(samFile, details, threads=1):
    #convert format to bam and index
    LOG.write("convertSamToBam(%s, details, %d)\n"%(samFile, threads))
    tempTime = time()
    sortThreads = max(int(threads/2), 1)
    bamFileUnsorted = re.sub(".sam", "_unsorted.bam", samFile, re.IGNORECASE)
    command = "samtools view -bS -@ %d -o %s %s"%(sortThreads, bamFileUnsorted, samFile)
    LOG.write("executing:\n"+command+"\n")
    return_code = subprocess.call(command, shell=True, stderr=LOG)
    LOG.write("samtools view return code = %d, time=%d\n"%(return_code, time()-tempTime))
    if return_code != 0:
        comment = "samtools view returned %d"%return_code
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None

    #os.remove(samFile) #save a little space
    bamFileSorted = re.sub(".sam", ".bam", samFile, re.IGNORECASE)
    command = "samtools sort -@ %d -o %s %s"%(sortThreads, bamFileSorted, bamFileUnsorted)
    LOG.write("executing:\n"+command+"\n")
    LOG.flush()
    return_code = subprocess.call(command, shell=True, stderr=LOG)
    LOG.write("samtools sort return code = %d, time=%d\n"%(return_code, time()-tempTime))
    if return_code != 0:
        comment = "samtools sort returned %d"%return_code
        LOG.write(comment+"\n")
        details["problem"].append(comment)
        return None

    command = ["samtools", "index", bamFileSorted]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    LOG.write("samtools index return code = %d\n"%return_code)
    return bamFileSorted

def runRacon(contigFile, longReadsFastq, details, threads=1):
    """
    Polish (correct) sequence of assembled contigs by comparing to the original long-read sequences
    Run racon on reads, read-to-contig-sam, contigs. Generate polished contigs.
    Return name of polished contigs.
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    LOG.write('runRacon(%s, %s, details, %d)\n'%(contigFile, longReadsFastq, threads))
    readsToContigsSam = runMinimap(contigFile, longReadsFastq, details, threads, outformat='sam')
    raconStartTime = time()
    raconContigs = contigFile.replace(".fasta", ".racon.fasta")
    raconOut = open(raconContigs, 'w')
    command = ["racon", "-t", str(threads), "-u", longReadsFastq, readsToContigsSam, contigFile]
    tempTime = time()
    LOG.write("racon command: \n"+' '.join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stderr=FNULL, stdout=raconOut)
    LOG.write("racon return code = %d, time = %d seconds\n"%(return_code, time()-raconStartTime))
    if return_code != 0:
        return None
    raconContigSize = os.path.getsize(raconContigs)
    if raconContigSize < 10:
        return None
    comment = "racon, input %s, output %s"%(contigFile, raconContigs)
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return raconContigs

def runBowtie(contigFile, shortReadFastq, details, threads=1, outformat='bam'):
    """
    index contigsFile, then run bowtie2, then convert sam file to pos-sorted bam and index
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    command = ["bowtie2-build", "--threads", str(threads), contigFile, contigFile]
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout and stderr to dev/null
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("bowtie2-build return code = %d\n"%return_code)
    if return_code != 0:
        return None

    command = ["bowtie2", "-p", str(threads), "-x", contigFile]
    fastqBase=''
    if ":" in shortReadFastq:
        read1, read2 = shortReadFastq.split(":")
        command.extend(('-1', read1, '-2', read2))
        fastqBase = read1
    else:
        command.extend(('-U', shortReadFastq))
        fastqBase = os.path.basename(shortReadFastq)
    fastqBase = re.sub(r"\..*", "", fastqBase)
    samFile = contigFile+"_"+fastqBase+".sam"
    command.extend(('-S', samFile))
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("bowtie2 return code = %d\n"%return_code)
    if return_code != 0:
        return None
    if outformat == 'sam':
        return samFile

    else:
        contigsBam = convertSamToBam(samFile, details, threads=threads)
        LOG.write('runBowtie returning %s\n'%contigsBam)
        return contigsBam

def runPilon(contigFile, shortReadFastq, details, pilon_jar, threads=1):
    """ 
    polish contigs with short reads (illumina or iontorrent)
    first map reads to contigs with bowtie
    """
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    if not os.path.exists(pilon_jar):
        comment = "jarfile %s not found when processing %s, giving up"%(pilon_jar, shortReadFastq)
        details['problem'].append(comment)
        LOG.write(comment+"\n")
        return

    bamFile = runBowtie(contigFile, shortReadFastq, details, threads=threads, outformat='bam')
    if not bamFile:
        return None
    command = ['java', '-Xmx32G', '-jar', pilon_jar, '--genome', contigFile]
    if ':' in shortReadFastq:
        command.extend(('--frags', bamFile))
    else:
        command.extend(('--unpaired', bamFile))
    pilonPrefix = contigFile.replace(".fasta", ".pilon")
    command.extend(('--outdir', '.', '--output', pilonPrefix, '--changes'))
    command.extend(('--threads', str(threads)))
    tempTime = time()
    LOG.write("executing:\n"+" ".join(command)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("pilon return code = %d\n"%return_code)
    LOG.write("pilon duration = %d\n"%(time() - tempTime))
    if return_code != 0:
        return None
    pilonContigs = pilonPrefix+".fasta"
    details['pilon_changes'] = 0
    with open(pilonContigs.replace(".fasta", ".changes")) as CHANGES:
        details['pilon_changes'] = len(CHANGES.read().splitlines())

    comment = "pilon, input %s, output %s, num_changes = %d"%(contigFile, pilonContigs, details['pilon_changes'])
    LOG.write(comment+"\n")
    details["post-assembly transformation"].append(comment)
    return pilonContigs 

def calcReadDepth(bamfiles):
    """ Return dict of contig_ids to tuple of (coverage, normalized_coverage) """
    LOG.write("calcReadDepth(%s)\n"%" ".join(bamfiles))
    readDepth = {}
    command = ["samtools", "depth", "-a"]
    if type(bamfiles) is str:
        command.append(bamfiles)
    else:
        command.extend(bamfiles)
    LOG.write("command = "+" ".join(command)+"\n")
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    depthData = proc.communicate()[0]
    LOG.write("length of depthData string = %d\n"%len(depthData))
    depthSum = 0
    totalDepthSum = 0
    totalLength = 0
    length = 0
    contigLength = {}
    prevContig=None
    for line in iter(depthData.splitlines()):
        fields = line.rstrip().split("\t")
        if len(fields) < 3:
            raise Exception("Number of fields is less than 3:\n"+line)
        contig = fields[0]
        depth = 0
        for field in fields[2:]:
            depth += float(field)
        if contig != prevContig:
            if prevContig is not None:
                meanDepth = depthSum/length
                readDepth[prevContig] = [meanDepth, 0]
                contigLength[prevContig] = length
                depthSum = 0
                length = 0
            prevContig = contig
        totalLength += 1
        length += 1
        depthSum += depth
        totalDepthSum += depth
    # process data for last contig
    if prevContig is not None:
        meanDepth = depthSum/length
        readDepth[prevContig] = [meanDepth, 0]
        contigLength[prevContig] = length

    LOG.write("total length for depth data = %d\n"%totalLength)
    LOG.write("total depth = %.1f\n"%totalDepthSum)
    LOG.write("len(readDepth) = %d\n"%len(readDepth))
    totalMeanDepth = 0
    if totalLength > 0:
        totalMeanDepth = totalDepthSum/totalLength

    # calculate mean depth of contigs within "normal" boundary around overall mean
    lowerBound = totalMeanDepth * 0.5
    upperBound = totalMeanDepth * 1.5
    oneXSum = 0.0
    oneXLen = 0
    for c in readDepth:
        meanDepth = readDepth[c][0]
        if meanDepth >= lowerBound and meanDepth <= upperBound:
            oneXSum += meanDepth * contigLength[c]
            oneXLen += contigLength[c]
    oneXDepth = 1
    if oneXLen > 0 and oneXSum > 0:
        oneXDepth = oneXSum/oneXLen # length-weighted average

    for c in readDepth:
        meanDepth = readDepth[c][0]
        normalizedDepth = meanDepth / oneXDepth
        readDepth[c][1] = normalizedDepth
    return readDepth

def runCanu(details, threads=1, genome_size="5m", memory=250, prefix=""):
    LOG.write("Time = %s, total elapsed = %d seconds\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time())), time()-START_TIME))
    canuStartTime = time()
    LOG.write("runCanu: elapsed seconds = %d\n"%(canuStartTime-START_TIME))
    comment = """
usage: canu [-version] [-citation] \
            [-correct | -trim | -assemble | -trim-assemble] \
            [-s <assembly-specifications-file>] \
            -p <assembly-prefix> \
            -d <assembly-directory> \
            genomeSize=<number>[g|m|k] \
            [other-options] \
            [-pacbio-raw | -pacbio-corrected | -nanopore-raw | -nanopore-corrected] file1 file2 ...
    """
    # canu -d /localscratch/allan/canu_assembly -p p6_25X gnuplotTested=true genomeSize=5m useGrid=false -pacbio-raw pacbio_p6_25X.fastq
    command = ["canu", "-d", '.', "-p", "canu", "useGrid=false", "genomeSize=%s"%genome_size]
    command.extend(["maxMemory=" + str(memory), "maxThreads=" + str(threads)])
    command.append("stopOnReadQuality=false")
    """
    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """
    pacbio_reads = []
    for item in details['reads']:
        if details['reads'][item]['platform'] in ('pacbio', 'fasta'):
            pacbio_reads.append(item)
            if details['reads'][item]['platform'] == 'fasta':
                comment = 'submitting fasta reads to canu, but calling them "pacbio": '+' '.join(details['platform']['fasta'])
                LOG.write(comment+"\n")
                details['problem'].append(comment)
    if pacbio_reads:
        command.append("-pacbio-raw")
        command.extend(pacbio_reads)
    nanopore_reads = []
    for item in details['reads']:
        if details['reads'][item]['platform'] == 'nanopore':
            nanopore_reads.append(item)
    if nanopore_reads:
        command.append("-nanopore-raw")
        command.extend(nanopore_reads)
    if not pacbio_reads + nanopore_reads:
        LOG.write("no long read files available for canu.\n")
        details["problem"].append("no long read files available for canu")
        return None
    LOG.write("canu command =\n"+" ".join(command)+"\n")

    canuStartTime = time()
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big
        return_code = subprocess.call(command, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("return code = %d\n"%return_code)
    canuEndTime = time()
    elapsedTime = canuEndTime - canuStartTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)

    details["assembly"] = { 
                'assembly_elapsed_time' : elapsedHumanReadable,
                'assembler': 'canu',
                'command_line': " ".join(command)
                }

    LOG.write("Duration of canu run was %s\n"%(elapsedHumanReadable))
    if os.path.exists("canu.report"):
        LOG.write("details_dir = %s\n"%DETAILS_DIR)
        LOG.write("canu_report file name = %s\n"%(prefix+"canu_report.txt"))
        canuReportFile = os.path.join(DETAILS_DIR, (prefix+"canu_report.txt"))
        LOG.write("moving canu.report to %s\n"%canuReportFile)
        shutil.move("canu.report", canuReportFile)
    
    if not os.path.exists("canu.contigs.fasta"):
        comment = "Canu failed to generate contigs file. Check "+prefix+"canu_report.txt"
        LOG.write(comment+"\n")
        details["assembly"]["outcome"] = comment
        details["problem"].append(comment)
        return None
    # rename to canonical contigs.fasta
    contigsFile = "contigs.fasta"
    shutil.move("canu.contigs.fasta", contigsFile)
    shutil.move("canu.contigs.gfa", os.path.join(DETAILS_DIR, prefix+"assembly_graph.gfa"))
    details["assembly"]["contigs.fasta size:"] = os.path.getsize(contigsFile)
    return contigsFile

def write_html_report(htmlFile, details):
    LOG.write("writing html report to %s\n"%htmlFile)
    HTML = open(htmlFile, 'w')
    HTML.write("<head><style>\n.a { left-margin: 50px }\n")
    HTML.write(".b {left-margin: 75px }\n")
    HTML.write("</style></head>\n")
    HTML.write("<h1>Genome Assembly Report</h1>\n")
    HTML.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")

    HTML.write("<h3>Input reads:</h3>\n")
    for item in details['reads']:
        if 'supercedes' in details['reads'][item]:
            continue # this is a derived item, not original input
        HTML.write(item+"<table class='a'>")
        for key in sorted(details['reads'][item]):
            if key == 'problem':
                continue
            HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['reads'][item][key])))
        HTML.write("</table>\n")
        if "problem" in details['reads'][item] and details['reads'][item]['problem']:
            HTML.write("<div class='b'><b>Issues with read set "+item+"</b>\n<ul>")
            for prob in details['reads'][item]['problem']:
                HTML.write("<li>"+prob+"\n")
            HTML.write("</ul></div>\n")
    
    if details["pre-assembly transformation"]:
        HTML.write("<h3>Pre-Assembly Transformations</h3>\n<div class='a'>\n")
        HTML.write("<ul>\n")
        for line in details["pre-assembly transformation"]:
            HTML.write("<li>"+line+"\n")
        HTML.write("</ul></div>\n")

    if "trim report" in details:
        HTML.write("<h3>Trimming Report</h3>\n<div class='a'>\n")
        for reads in details["trim report"]:
            HTML.write("<b>"+reads+"</b><ul>")
            for report in details["trim report"][reads]:
                if os.path.exists(os.path.join(DETAILS_DIR, report)):
                    HTML.write("<pre>\n")
                    HTML.write(open(os.path.join(DETAILS_DIR, report)).read())
                    HTML.write("\n</pre>\n")
        HTML.write("</div>\n")

    if details['derived_reads']:
        HTML.write("<h3>Transformed reads:</h3>\n")
        for item in details['derived_reads']:
            HTML.write(item+"<table class='a'>")
            for key in sorted(details['reads'][item]):
                if key == 'problem':
                    continue
                HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['reads'][item][key])))
            HTML.write("</table>\n")
            if "problem" in details['reads'][item] and details['reads'][item]['problem']:
                HTML.write("<div class='b'><b>Issues with read set "+item+"</b>\n<ul>")
                for prob in details['reads'][item]['problem']:
                    HTML.write("<li>"+prob+"\n")
                HTML.write("</ul></div>\n")

    HTML.write("<h3>Assembly</h3>\n")
    if 'assembly' in details:
        HTML.write("<table class='a'>")
        for key in sorted(details['assembly']):
            HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(details['assembly'][key])))
        HTML.write("</table>\n")
    else:
        HTML.write("<p>None</p>\n")

    if "quast_txt" in details:
        HTML.write("<h3>Quast Report</h3>\n")
        HTML.write("<table class='a'>")
        HTML.write("<li><a href='%s'>%s</a>\n"%(details["quast_txt"], "Quast text report"))
        HTML.write("<li><a href='%s'>%s</a>\n"%(details["quast_html"], "Quast html report"))
        HTML.write("</table>\n")
        if os.path.exists(os.path.join(SAVE_DIR, details["quast_txt"])):
            HTML.write("<pre>\n")
            HTML.write(open(os.path.join(SAVE_DIR, details["quast_txt"])).read())
            HTML.write("\n</pre>\n")
    
    if details["post-assembly transformation"]:
        HTML.write("<h3>Post-Assembly Transformations</h3>\n<div class='a'>\n<ul>\n")
        HTML.write("<ul>\n")
        for line in details['post-assembly transformation']:
            HTML.write("<li>"+line+"\n")
        HTML.write("</ul></div>\n")

    if "circular_contigs" in details and details['circular_contigs']:
        HTML.write("<h3>Circular Contigs</h3>\n<div class='a'>\n")
        HTML.write("<b>As suggested by Unicycler</b>\n")
        HTML.write("<ul>\n")
        for line in details['circular_contigs']:
            HTML.write("<li>"+line+"\n")
        HTML.write("</ul></div>\n")

    if "contig_filtering" in details and details['contig_filtering']:
        HTML.write("<h3>Contig Filtering</h3>\n<div class='a'>\n")
        HTML.write("<ul>\n")
        for line in details['contig_filtering']:
            HTML.write("<li>"+line+"\n")
        HTML.write("</ul>\n")

        if "bad_contigs" in details and details["bad_contigs"]:
            HTML.write("<b>Contigs Below Thresholds</b>\n")
            if len(details["bad_contigs"]) > 10:
                HTML.write("<b>Showing first 10 entries</b>\n")
            HTML.write("<ul>\n")
            for line in details['bad_contigs'][:10]:
                HTML.write("<li>"+line+"\n")
            HTML.write("</ul>\n")
        HTML.write("</div>\n")

    if "Bandage plot" in details:
        #path, imageFile = os.path.split(details["Bandage plot"])
        if os.path.exists(details["Bandage plot"]):
            svg_text = open(details["Bandage plot"]).read()
            svg_text = re.sub(r'<svg width="[\d\.]+mm" height="[\d\.]+mm"', '<svg width="200mm" height="150mm"', svg_text)
            HTML.write("<h3>Bandage Plot</h3>\n")
            HTML.write("<div class='a'>")
            HTML.write(svg_text+"\n\n")
            #HTML.write("<img src='%s'>\n"%imageFile)
            HTML.write("</div>\n")
    HTML.write("</html>\n")
    HTML.close()


def main():
    global START_TIME
    START_TIME = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='p3_assembly')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', metavar='files', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', metavar='files', nargs='*', help='list of IonTorrent[.gz] files or pairs, use : between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', metavar='files', nargs='*', help='list of Pacific Biosciences fastq[.gz] files', required=False, default=[])
    parser.add_argument('--nanopore', metavar='files', nargs='*', help='list of Oxford Nanotech fastq[.gz] files', required=False, default=[])
    parser.add_argument('--sra', metavar='files', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', metavar='files', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('--recipe', choices=['unicycler', 'canu', 'spades', 'meta-spades', 'plasmid-spades', 'single-cell', 'auto'], help='assembler to use', default='auto')

    parser.add_argument('--racon_iterations', type=int, default=2, help='number of times to run racon per long-read file', required=False)
    parser.add_argument('--pilon_iterations', type=int, default=2, help='number of times to run pilon per short-read file', required=False)
    #parser.add_argument('--singlecell', action = 'store_true', help='flag for single-cell MDA data for SPAdes', required=False)
    parser.add_argument('--prefix', default='', help='prefix for output files', required=False)
    parser.add_argument('--genome_size', metavar='k, m, or g', default=DEFAULT_GENOME_SIZE, help='genome size for canu: e.g. 300k or 5m or 1.1g', required=False)
    parser.add_argument('--min_contig_length', type=int, default=300, help='save contigs of this length or longer', required=False)
    parser.add_argument('--min_contig_coverage', type=float, default=5, help='save contigs of this coverage or deeper', required=False)
    parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--no_pilon', action='store_true', help='for unicycler', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('-t', '--threads', metavar='cores', type=int, default=4)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, default=250, help='RAM limit in Gb')
    parser.add_argument('--trim', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--pilon_jar', help='path to pilon executable or jar')
    parser.add_argument('--bandage', action='store_true', help='generate image of assembly path using Bandage')
    parser.add_argument('--params_json', help='JSON file with additional information.')
    parser.add_argument('--path-prefix', help="Add the given directories to the PATH", nargs='*', required=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    if args.params_json:
        parseJsonParameters(args)
    baseName = args.outputDirectory #"p3_assembly" 
    if args.prefix and not args.prefix.endswith("_"):
        args.prefix += "_"
    global WORK_DIR
    WORK_DIR = baseName+"_work"
    if os.path.exists(WORK_DIR):
        shutil.rmtree(WORK_DIR)
    os.mkdir(WORK_DIR)
    WORK_DIR = os.path.abspath(WORK_DIR)
    global SAVE_DIR
    SAVE_DIR = os.path.abspath(os.path.join(WORK_DIR, "save"))
    if os.path.exists(SAVE_DIR):
        shutil.rmtree(SAVE_DIR)
    os.mkdir(SAVE_DIR)
    global DETAILS_DIR
    DETAILS_DIR = os.path.abspath(os.path.join(SAVE_DIR, "details"))
    os.mkdir(DETAILS_DIR)
    logfileName = os.path.join(DETAILS_DIR, args.prefix + "p3_assembly.log")
 
    if args.path_prefix:
        os.environ["PATH"] = ":".join(args.path_prefix) + ":" + os.environ["PATH"]

    global LOG 
    sys.stderr.write("logging to "+logfileName+"\n")
    LOG = open(logfileName, 'w', 0) #unbuffered 
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(START_TIME))+"\n")
    LOG.write("args= "+str(args)+"\n\n")
    LOG.write("Work directory is "+WORK_DIR+"\n\n")
    LOG.write("Final output will be saved to "+SAVE_DIR+"\n\n")
    LOG.write("Detailed output will be saved to "+DETAILS_DIR+"\n\n")
    details = { 'logfile' : logfileName }
    details["pre-assembly transformation"] = []
    details["post-assembly transformation"] = []
    details["original_items"] = []
    details["reads"] = {}
    details["problem"] = []
    details["derived_reads"] = []
    details["platform"] = {'illumina':[], 'iontorrent':[], 'pacbio':[], 'nanopore':[], 'fasta':[], 'anonymous':[]}

    if args.illumina:
        platform = 'illumina'
        for item in args.illumina:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(item, details, platform=platform, interleaved=interleaved)

    if args.iontorrent:
        platform = 'iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            registerReads(item, details, platform=platform, interleaved=interleaved)

    if args.pacbio:
        for item in args.pacbio:
            registerReads(item, details, platform = 'pacbio')

    if args.nanopore:
        for item in args.nanopore:
            registerReads(item, details, platform = 'nanopore')

    if args.fasta:
        for item in args.fasta:
            registerReads(item, details, platform = 'fasta')

    if args.sra:
        fetch_sra_files(args.sra, details)

    if args.anonymous_reads:
        categorize_anonymous_read_files(args, details)

    # move into working directory so that all files are local
    os.chdir(WORK_DIR)

    if args.trim and details['platform']['illumina'] + details['platform']['iontorrent']:
        trimGalore(details, threads=args.threads)
    LOG.write("details dir = "+DETAILS_DIR+"\n")
    if args.recipe == "auto":
        #now must decide which assembler to use
        if True:
            # original rule: if any illumina or iontorrent reads present, use Unicycler (long-reads can be present), else use canu for long-reads
            if details['platform']['illumina'] + details['platform']['iontorrent']:
                args.recipe = "unicycler"
            else:
                args.recipe = "canu"
        else:
            # alternative rule: if any long reads present, use canu
            if details['platform']['pacbio'] + details['platform']['nanopore']:
                args.recipe = "canu"
            else:
                args.recipe = "unicycler"
    if "spades" in args.recipe or args.recipe == "single-cell":
        contigs = runSpades(details, args)
    elif args.recipe == "unicycler":
        contigs = runUnicycler(details, threads=args.threads, min_contig_length=args.min_contig_length, prefix=args.prefix)
    elif args.recipe == "canu":
        contigs = runCanu(details, threads=args.threads, genome_size=args.genome_size, memory=args.memory, prefix=args.prefix)
    else:
        LOG.write("cannot interpret args.recipe: "+args.recipe)

    if contigs and os.path.getsize(contigs):
        # now run racon with each long-read file
        for longReadFile in details['reads']:
            if details['reads'][longReadFile]['length_class'] == 'long':
                for i in range(0, args.racon_iterations):
                    LOG.write("runRacon(%s, %s, details, threads=%d)\n"%(contigs, longReadFile, args.threads))
                    raconContigFile = runRacon(contigs, longReadFile, details, threads=args.threads)
                    if raconContigFile is not None:
                        contigs = raconContigFile
                    else:
                        break # break out of iterating racon_iterations, go to next long-read file if any
        
    if contigs and os.path.getsize(contigs):
        # now run pilon with each short-read file
        for shortReadFastq in details['reads']:
            if 'superceded_by' in details['reads'][shortReadFastq]:
                continue # may have been superceded by trimmed version of those reads
            if details['reads'][shortReadFastq]['length_class'] == 'short':
                for iteration in range(0, args.pilon_iterations):
                    LOG.write("runPilon(%s, %s, details, %s, threads=%d) iteration=%d\n"%(contigs, shortReadFastq, args.pilon_jar, args.threads, iteration))
                    pilonContigFile = runPilon(contigs, shortReadFastq, details, args.pilon_jar, threads=args.threads)
                    if pilonContigFile is not None:
                        contigs = pilonContigFile
                    else:
                        break
                    if details['pilon_changes'] == 0:
                        break
        
    if contigs and os.path.getsize(contigs):
        filteredContigs = filterContigsByMinLength(contigs, details, args.min_contig_length, args.min_contig_coverage, args.threads, args.prefix)
        if filteredContigs:
            contigs = filteredContigs
    if contigs and os.path.getsize(contigs):
        runQuast(contigs, args, details)
        shutil.move(contigs, os.path.join(SAVE_DIR, args.prefix+"contigs.fasta"))

    gfaFile = os.path.join(DETAILS_DIR, args.prefix+"assembly_graph.gfa")
    if os.path.exists(gfaFile):
        bandagePlot = runBandage(gfaFile, details)
        details["Bandage plot"] = bandagePlot

    with open(os.path.join(DETAILS_DIR, args.prefix+"run_details.json"), "w") as fp:
        json.dump(details, fp, indent=2, sort_keys=True)
    htmlFile = os.path.join(SAVE_DIR, args.prefix+"assembly_report.html")
    write_html_report(htmlFile, details)
    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    LOG.write("Total time in hours = %d\t"%((time() - START_TIME)/3600))
    LOG.close()


if __name__ == "__main__":
    main()
