#!/usr/bin/env python
import sys
import subprocess
import argparse
import os
import os.path
import re
import shutil
from time import time, localtime, strftime
import json
import glob
from ReadLibrary import ReadLibrary
from ReadLibrary import constrain_total_bases

"""
This script organizes a command line for an assembly program: 
Flye, Unicycler, Spades or canu as appropriate: 
    canu if only long reads (pacbio or nanopore), 
    Unicycler if illumina or iontorrent 
    and Unicycler for hybrid assemblies, eg Illumina plus PacBio
    or Spades if requested
It can auto-detect read types (illumina, iontorrent, pacbio, nanopore)
It can run trim_galore prior to assembling.
It can run Quast to generate assembly quality statistics.
It can perform 'polishing' using Pilon (for Illumina) or Racon (for pacbio or nanopore reads)
TODO: properly handle different kinds of pacbio reads
"""

DEFAULT_GENOME_SIZE = 5000000
DEFAULT_MAX_BASES=1e10 # this is 2000X coverage for a 5M genome
DEFAULT_MEMORY_GIGABYTES = 128
LOG = None # create a log file at start of main()
START_TIME = None
WORK_DIR = None
SAVE_DIR = None
DETAILS_DIR = None

def runQuast(contigsFile, min_contig_length, num_threads=1):
    LOG.write("runQuast: time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    quastDir = "quast_out"
    quastCommand = ["quast.py",
                    "-o", quastDir,
                    "-t", str(num_threads),
                    "--min-contig", str(min_contig_length),
                    contigsFile]
    LOG.write("running quast: "+" ".join(quastCommand)+"\n")
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null
        return_code = subprocess.call(quastCommand, shell=False, stdout=FNULL, stderr=FNULL)
    LOG.write("return code = %d\n"%return_code)
    report = {}
    if return_code == 0:
        report['program'] = 'quast'
        report['command'] = ' '.join(quastCommand)
        report['report_file_path'] = os.path.join(quastDir, "report.txt")
        with open(os.path.join(quastDir, "report.txt")) as F:
            report['report_text'] = F.read()
        quastCommand = ["quast.py", "--version"]
        proc = subprocess.Popen(quastCommand, shell=False, stdout=subprocess.PIPE)
        version_text = proc.stdout.read().decode()
        proc.wait()
        report["version"] = version_text
        return report
    else:
        return None

def filterContigsByLengthAndCoverage(inputContigs, read_list, min_contig_length=300, min_contig_coverage=5, num_threads=1, prefix=""):
    """ 
    Write only sequences at or above min_length and min coverage to output file.
    """
    LOG.write("filterContigsByLengthAndCoverage: Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    report = {}
    shortReadDepth = None
    longReadDepth = None
    bamFiles = []
    for read_set in read_list:
        if read_set.length_class == 'short':
            bowtieReport = runBowtie(inputContigs, read_set, num_threads, outformat='bam')
            if 'bam' in bowtieReport:
                bamFiles.append(bowtieReport['bam'])
                report['bowtie2_version'] = bowtieReport['version']
    if bamFiles:
        (average_depth, shortReadDepth) = calcReadDepth(bamFiles)
        report['average depth (short reads)'] = "{:.2f}".format(average_depth)
        # there are two calculations of read coverage depth, redundant but slightly different
    bamFiles = []
    for read_set in read_list:
        if read_set.length_class == 'long':
            minimapReport = runMinimap(inputContigs, read_set, num_threads, outformat='bam')
            if 'bam' in minimapReport:
                bamFiles.append(minimapReport['bam'])
                report['minimap_version'] = minimapReport['version']
    if bamFiles:
        (average_depth, longReadDepth) = calcReadDepth(bamFiles)
        report['average depth (long reads)'] = "{:.2f}".format(average_depth)
        # there are two calculations of read coverage depth, redundant but slightly different
    report['min_contig_length_threshold'] = "%d"%min_contig_length
    report["min_contig_coverage_threshold"] = "%.1f"%min_contig_coverage
    num_good_contigs = num_bad_contigs = 0
    suboptimalContigsFile = "contigs_below_length_coverage_threshold.fasta"
    total_seq_length = 0
    weighted_short_read_coverage = 0
    weighted_long_read_coverage = 0
    outputContigs = re.sub(r"\..*", "_depth_cov_filtered.fasta", inputContigs)
    LOG.write("writing filtered contigs to %s\n"%outputContigs)
    with open(inputContigs) as IN:
        with open(outputContigs, 'w') as OUT:
            SUBOPT = open(os.path.join(DETAILS_DIR, suboptimalContigsFile), "w")
            seqId=None
            seq = ""
            contigIndex = 1
            num_circular_contigs = 0
            line = "1"
            while line:
                line = IN.readline()
                m = re.match(r">(\S+)", line)
                if m or not line: 
                    if seq:
                        contigId = ">"+prefix+"contig_%d"%contigIndex
                        contigInfo = " length %5d"%len(seq)
                        contigIndex += 1
                        short_read_coverage = 0
                        long_read_coverage = 0
                        passes_thresholds = False
                        if shortReadDepth and seqId in shortReadDepth:
                            short_read_coverage, normalizedDepth = shortReadDepth[seqId]
                            contigInfo += " coverage %.01f normalized_cov %.2f"%(short_read_coverage, normalizedDepth)
                            passes_thresholds = short_read_coverage >= min_contig_coverage
                        if longReadDepth and seqId in longReadDepth:
                            long_read_coverage, normalizedDepth = longReadDepth[seqId]
                            contigInfo += " longread_coverage %.01f normalized_longread_cov %.2f"%(long_read_coverage, normalizedDepth)
                            passes_thresholds |= long_read_coverage >= min_contig_coverage
                        if len(seq) < min_contig_length:
                            passes_thresholds = False
                            num_bad_contigs += 1
                        if passes_thresholds:
                            OUT.write(contigId+contigInfo+"\n")
                            for i in range(0, len(seq), 60):
                                OUT.write(seq[i:i+60]+"\n")
                            num_good_contigs += 1
                            if short_read_coverage:
                                weighted_short_read_coverage += short_read_coverage * len(seq)
                            if long_read_coverage:
                                weighted_long_read_coverage += long_read_coverage * len(seq)
                            total_seq_length += len(seq)
                            if "circular=true" in line:
                                num_circular_contigs += 1
                        else:
                            SUBOPT.write(contigId+contigInfo+"\n")
                            for i in range(0, len(seq), 60):
                                SUBOPT.write(seq[i:i+60]+"\n")
                            num_bad_contigs += 1
                        seq = ""
                    if m:
                        seqId = m.group(1)
                elif line:
                    seq += line.rstrip()
    if False and total_seq_length:
        if weighted_short_read_coverage:
            report['average short read coverage'] = "%.3f"%(weighted_short_read_coverage / total_seq_length)
        if weighted_long_read_coverage:
            report['average long read coverage'] = "%.3f"%(weighted_long_read_coverage / total_seq_length)
    report['total contigs'] = "%d"%(num_good_contigs+num_bad_contigs)
    report['num contigs above thresholds'] = "%d"%num_good_contigs
    report['num contigs below thresholds'] = "%d"%num_bad_contigs
    report['total length of good contigs'] = "%d"%total_seq_length
    if num_circular_contigs:
        report['contigs predicted circular'] = "%d"%num_circular_contigs
    if num_bad_contigs:
        report["suboptimal contigs file"] = suboptimalContigsFile
    if num_good_contigs:
        report['good contigs file'] = outputContigs
    else:
        LOG.write("failed to generate outputContigs\n")
    comment = "filterContigsByLengthAndCoverage, input %s, output %s"%(inputContigs, outputContigs)
    LOG.write(comment+"\n")
    return report

def runBandage(gfaFile):
    imageFormat = ".svg"
    report = {}
    if os.path.exists(gfaFile):
        plotFile = gfaFile.replace(".gfa", ".plot"+imageFormat)
        command = ["Bandage", "image", gfaFile, plotFile]
        LOG.write(" ".join(command)+"\n")
        report['command'] = ' '.join(command)
        try:
            with open(os.devnull, 'w') as FNULL:
                return_code = subprocess.call(command, shell=False, stderr=FNULL)
            LOG.write("return code = %d\n"%return_code)
            if return_code == 0:
                retval = plotFile
                proc = subprocess.Popen(["Bandage", "--version"], shell=False, stdout=subprocess.PIPE)
                version_text = proc.stdout.read().decode().strip()
                proc.wait()
                report['program'] = 'Bandage'
                report['version'] = version_text
                report['plot'] = plotFile
            else:
                LOG.write("Error creating Bandage plot\n")
        except OSError as ose:
            comment = "Problem running Bandage: "+str(ose)
            LOG.write(comment+"\n")
    return report

def runMegahit(read_list, prefix, threads=4, memoryGB=16):
    LOG.write("runMegahit\n")
    proc = subprocess.Popen(["megahit", "--version"], shell=False, stdout=subprocess.PIPE)
    version_text = proc.stdout.read().decode()
    if version_text:
        m = re.search(r"MEGAHIT\s+\S+", version_text, flags=re.IGNORECASE)
        if m:
            version_text = m.group(0) # entire match
    proc.wait()
    assembly_data = {}
    assembly_data['program'] = 'megahit'
    assembly_data['version'] = version_text

    memoryBytes = memoryGB * 1024 * 1024 *1024
    command = ["megahit", "-t", str(threads), '-m', str(memoryBytes)]
    for read_set in read_list:
        if read_set.length_class == 'short':
            if len(read_set.files) > 1:
                command.extend(("-1", read_set.files[0], "-2", read_set.files[1]))
            else:
                command.extend(("--read", read_set.files[0]))

    LOG.write(" ".join(command)+"\n")
    LOG.flush()
    assembly_data['command'] = " ".join(command)
    proc = subprocess.run(command, shell=False)
    proc.check_returncode()
    #LOG.write("return code = %d\n"%return_code)

    if os.path.exists("megahit_out/log"):
        megahitLogFile = prefix+"megahit.log"
        shutil.move("megahit_out/log", os.path.join(DETAILS_DIR, megahitLogFile))

    if os.path.exists("megahit_out/final.contigs.fa"):
        shutil.move("megahit_out/final.contigs.fa", "contigs.fasta") #rename to canonical name
    return assembly_data

def runUnicycler(read_list, threads=1, min_contig_length=0, prefix="", spades_exec=None):
    LOG.write("runUnicycler: Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    proc = subprocess.Popen(["unicycler", "--version"], shell=False, stdout=subprocess.PIPE)
    version_text = proc.stdout.read().decode()
    if version_text:
        m = re.search(r"Unicycler\s+\S+", version_text, flags=re.IGNORECASE)
        if m:
            version_text = m.group(0) # entire match
    proc.wait()
    assembly_data = {}
    assembly_data['program'] = 'unicycler'
    assembly_data['version'] = version_text

    command = ["unicycler", "-t", str(threads), "-o", '.']
    if min_contig_length:
        command.extend(("--min_fasta_length", str(min_contig_length)))
    command.extend(("--keep", "2")) # keep files needed for re-run if necessary
    command.append("--no_pilon")  # we will run our own, if requested
    if spades_exec:
        command.extend(("--spades_path", spades_exec));

    # apparently unicycler can only accept one read set in each class (I tried multiple ways to submit 2 paired-end sets, failed)
    for read_set in read_list:
        if read_set.length_class == 'short':
            if len(read_set.files) > 1:
                command.extend(("--short1", read_set.files[0], "--short2", read_set.files[1]))
            else:
                command.extend(("--unpaired", read_set.files[0]))

        else:
            command.extend(("--long", read_set.files[0]))

    LOG.write(" ".join(command)+"\n")
    LOG.flush()
    with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big and unicycle.log is better
        return_code = subprocess.call(command, shell=False, stdout=FNULL)
    LOG.write("return code = %d\n"%return_code)

    if not (os.path.exists("assembly.fasta") and os.path.getsize("assembly.fasta")):
        comment = "First run of Unicycler resulted in no assembly, try again with more lenient parameters."
        LOG.write(comment+"\n")
        command.extend(("--mode", "bold", "--min_component_size", "300", "--min_dead_end_size", "300", "--depth_filter", "0.1"))
        comment = "re-run unicycler with command = "+" ".join(command)
        LOG.write(comment+"\n")
        with open(os.devnull, 'w') as FNULL: # send stdout to dev/null, it is too big and unicycle.log is better
            return_code = subprocess.call(command, shell=False, stdout=FNULL)
        LOG.write("return code = %d\n"%return_code)

    assembly_data['command'] = " ".join(command)
    unicyclerLogFile = "unicycler.log"
    if os.path.exists("unicycler.log"):
        unicyclerLogFile = prefix+"unicycler.log"
        shutil.move("unicycler.log", os.path.join(DETAILS_DIR, unicyclerLogFile))

    if not os.path.exists("assembly.fasta"):
        comment = "Unicycler failed to generate assembly file. Check "+unicyclerLogFile
        LOG.write(comment+"\n")
        return None

    assemblyGraphFile = prefix+"assembly_graph.gfa"
    shutil.move("assembly.gfa", os.path.join(DETAILS_DIR, assemblyGraphFile))

    contigsFile = "contigs.fasta"
    shutil.move("assembly.fasta", contigsFile) #rename to canonical name
    return assembly_data

def runSpades(read_list, prefix="", recipe=None, threads=4, memoryGB=DEFAULT_MEMORY_GIGABYTES):
    LOG.write("runSpades Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime())))
    command = ["spades.py", "--version"]
    proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    version_text = proc.stdout.read().decode().strip()
    if not version_text: # as of v3.14.1 this output came to stderr instead of stdout
        version_text = proc.stderr.read().decode().strip()
    proc.wait()
    assembly_data = {}
    assembly_data['program'] = 'spades'
    assembly_data['version'] = version_text

    command = ["spades.py", "--threads", str(threads), "-o", "."]
    if recipe == 'single-cell':
        command.append("--sc")
    if memoryGB:
        command.extend(["-m", str(memoryGB)])
    if recipe == "meta-spades":
        command.append("--meta")
    if recipe == "rna-spades":
        command.append("--rna")
        #
        # Validate arguments for metagenomic spades. It can only run with
        # a single paired-end library.
        #
        #if len(single_end_reads) > 0 or len(paired_end_reads[0]) > 1:
        #    sys.stderr.write("SPAdes in metagenomics mode can only process a single paired-end read file\n")
        #    sys.exit(1);
    if recipe == "plasmid-spades":
        command.append("--plasmid")
   
    any_fasta = False
    any_illumina = False
    any_iontorrent = False
    paired_end_counter = 1
    single_end_counter = 1
    for read_set in read_list:
        if read_set.format == "fasta":
            any_fasta = True
        if read_set.length_class == "short":
            if len(read_set.files) > 1:
                if paired_end_counter > 9:
                    LOG.write("Spades cannot take more than 9 paired-end libraries.")
                    continue 
                command.extend(("--pe{}-1".format(paired_end_counter), read_set.files[0], "--pe{}-2".format(paired_end_counter), read_set.files[1]))
                paired_end_counter += 1
            else:
                if single_end_counter > 9:
                    LOG.write("Spades cannot take more than 9 single-end libraries.")
                    continue 
                command.extend(("--s{}".format(single_end_counter), read_set.files[0]))
                single_end_counter += 1

        else: # length_class is long
            if read_set.platform == 'pacbio':
                command.extend(["--pacbio", read_set.files[0]])
            if read_set.platform == 'nanopore':
                command.extend(["--nanopore", read_set.files[0]])
            
    if any_fasta: # lacking quality scores means we need to turn off read correction
        command.append("--only-assembler") 
    if any_iontorrent:
        command.append("--iontorrent") # tell SPAdes that this is the read type
    LOG.write("SPAdes command =\n"+" ".join(command)+"\n")
    #LOG.write("    PATH:  "+os.environ["PATH"]+"\n\n")
    LOG.flush()

    return_code = subprocess.call(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    LOG.write("return code = %d\n"%return_code)

    contigsFile = "contigs.fasta"
    if return_code and not os.path.exists(contigsFile):
        comment = "spades failed to generate assembly, return code = %d"%return_code
        LOG.write(comment+"\n")
        comment = "try again adding '--only-assembler' option to spades"
        LOG.write(comment+"\n")
        command.append("--only-assembler")
        return_code = subprocess.call(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        LOG.write("return code = %d\n"%return_code)

    assembly_data['command'] = " ".join(command)
    spadesLogFile = "spades.log"
    if os.path.exists("spades.log"):
        shutil.move("spades.log", os.path.join(DETAILS_DIR, "spades.log"))
    if not os.path.exists(contigsFile):
        comment = "SPAdes failed to generate contigs file. Check "+spadesLogFile
        LOG.write(comment+"\n")
        return None
    assemblyGraphFile = prefix+"assembly_graph.gfa"
    if os.path.exists("assembly_graph_with_scaffolds.gfa"):
        LOG.write("found gfa file: assembly_graph_with_scaffolds.gfa\n")
        shutil.move("assembly_graph_with_scaffolds.gfa", os.path.join(DETAILS_DIR, assemblyGraphFile))
    else:
        gfa_candidates = glob.glob("*gfa")
        if gfa_candidates:
            LOG.write("found gfa file: {}\n".format(gfa_candidates[-1]))
            shutil.move(gfa_candidates[-1], os.path.join(DETAILS_DIR, assemblyGraphFile))
    return assembly_data

def runCanu(read_list, canu_exec="canu", threads=1, genome_size=DEFAULT_GENOME_SIZE, memory=DEFAULT_MEMORY_GIGABYTES, prefix=""):
    LOG.write("runCanu: Tiime = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    comment = """
usage: canu [-version] [-citation] \
            [-correct | -trim | -assemble | -trim-assemble] \
            [-s <assembly-specifications-file>] \
            -p <assembly-prefix> \
            -d <assembly-directory> \
            genomeSize=<number>[g|m|k] \
            [other-options] \
            [-pacbio-raw | -pacbio-corrected | -nanopore-raw | -nanopore-corrected] file1 file2 ...

    https://canu.readthedocs.io/en/latest/parameter-reference.html
    """
    # canu -d /localscratch/allan/canu_assembly -p p6_25X gnuplotTested=true genomeSize=5m useGrid=false -pacbio-raw pacbio_p6_25X.fastq
    # first get canu version
    p = subprocess.Popen([canu_exec, "--version"], shell=False, stdout=subprocess.PIPE)
    canu_version = p.stdout.readline().decode().rstrip()
    assembly_data = {}
    assembly_data['program'] = 'canu'
    assembly_data['version'] = canu_version
    p.wait()

    command = [canu_exec, "-d", '.', "-p", "canu", "genomeSize=%d"%genome_size]
    command.extend(["maxMemory=" + str(memory), "maxThreads=" + str(threads)])
    if "1.7" in canu_version:
        # special handling for this version
        command.append("gnuplotTested=true")
        command.append("stopOnReadQuality=false")
    
    pacbio_reads = []
    nanopore_reads = []
    for read_set in read_list:
        if read_set.platform == 'pacbio':
            pacbio_reads.append(read_set.files[0])
        if read_set.platform == 'nanopore':
            nanopore_reads.append(read_set.files[0])
    if pacbio_reads:
        command.append("-pacbio-raw")
        command.extend(pacbio_reads)
    if nanopore_reads:
        command.append("-nanopore-raw")
        command.extend(nanopore_reads)
    if not len(pacbio_reads) + len(nanopore_reads):
        LOG.write("no long read files available for canu.\n")
        return None
    LOG.write("canu command =\n"+" ".join(command)+"\n")
    LOG.flush()

    return_code = subprocess.call(command, shell=False)
    assembly_data['command'] = " ".join(command)

    if os.path.exists("canu.report"):
        LOG.write("details_dir = %s\n"%DETAILS_DIR)
        LOG.write("canu_report file name = %s\n"%(prefix+"canu_report.txt"))
        canuReportFile = os.path.join(DETAILS_DIR, (prefix+"canu_report.txt"))
        LOG.write("moving canu.report to %s\n"%canuReportFile)
        shutil.move("canu.report", canuReportFile)
    
    if not os.path.exists("canu.contigs.fasta"):
        comment = "Canu failed to generate contigs file. Check "+prefix+"canu_report.txt"
        LOG.write(comment+"\n")
        return None
    # rename to canonical contigs.fasta
    contigsFile = "contigs.fasta"
    shutil.move("canu.contigs.fasta", contigsFile)
    if os.path.exists("canu.contigs.gfa"):
        shutil.move("canu.contigs.gfa", os.path.join(DETAILS_DIR, prefix+"assembly_graph.gfa"))
    elif os.path.exists("canu.unitigs.gfa"):
        shutil.move("canu.unitigs.gfa", os.path.join(DETAILS_DIR, prefix+"assembly_graph.gfa"))
    return assembly_data

def runFlye(read_list, recipe="flye", threads=1, genome_size=DEFAULT_GENOME_SIZE, prefix="flye_", memory=128):
    LOG.write("runFlye: Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    comment = """
usage: flye (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw |
    --nano-corr | --nano-hq ) file1 [file_2 ...]
    --out-dir PATH

    [--genome-size SIZE] [--threads int] [--iterations int]
    [--meta] [--polish-target] [--min-overlap SIZE]
    [--keep-haplotypes] [--debug] [--version] [--help] 
    [--scaffold] [-uresume] [--resume-from] [--stop-after] 
    [--read-error float] [--extra-params]
    """
    # first get flye version
    p = subprocess.Popen(['flye', "--version"], shell=False, stdout=subprocess.PIPE)
    version_text = p.stdout.readline().decode().rstrip()
    assembly_data = {}
    assembly_data['program'] = 'flye'
    assembly_data['version'] = version_text
    p.wait()

    pacbio_reads = []
    nanopore_reads = []
    for read_set in read_list:
        #print("look at {}, type {}".platform(read_set.files[0], read_set.platform))
        if read_set.platform == 'pacbio':
            pacbio_reads.append(read_set.files[0])
        if read_set.platform == 'nanopore':
            nanopore_reads.append(read_set.files[0])
    command = ['flye', "--out-dir", '.', "--genome-size", str(genome_size), '--threads', str(threads)]
    if recipe == 'meta-flye':
        command.append("--meta")
    if pacbio_reads:
        command.append("--pacbio-raw")
        command.extend(pacbio_reads)
    elif nanopore_reads:
        command.append("--nano-raw")
        command.extend(nanopore_reads)
    #if memory:
    #    command.extend(['-str(memory)+"g")
    LOG.write(" ".join(command)+"\n")
    if not pacbio_reads + nanopore_reads:
        LOG.write("no long read files available for flye.\n")
        return None
    with open(os.path.join(DETAILS_DIR, prefix+"flye_stdout.txt"), 'w') as FLYE_STDOUT: 
        return_code = subprocess.call(command, shell=False, stdout=FLYE_STDOUT, stderr=FLYE_STDOUT)

    assembly_data['command'] = " ".join(command)

    if os.path.exists("flye.report"):
        LOG.write("flye_report file name = %s\n"%(prefix+"flye_report.txt"))
        flyeReportFile = os.path.join(DETAILS_DIR, (prefix+"flye_report.txt"))
        LOG.write("moving flye.report to %s\n"%flyeReportFile)
        shutil.move("flye.report", flyeReportFile)
    
    # rename to canonical contigs.fasta
    if os.path.exists("assembly.fasta"):
        shutil.move("assembly.fasta", "contigs.fasta")
    if os.path.exists("assembly_graph.gfa"):
        shutil.move("assembly_graph.gfa", os.path.join(DETAILS_DIR, prefix+"assembly_graph.gfa"))
    return assembly_data

def runMinimap(contigFile, read_set, num_threads, outformat='sam', save_intermediate=False):
    LOG.write("runMinimap: Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    """
    Map long reads to contigs by minimap2 (read paf-file, readsFile; generate paf file).
    """
    report = {}
    LOG.write("runMinimap(%s, %s, %s, %s)\n"%(contigFile, read_set.files[0], str(num_threads), outformat))
    command = ["minimap2", "--version"]
    proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, text=True)
    proc.wait()
    version_text = proc.stdout.read()
    report['program'] = 'minimap2'
    report['version'] = version_text.strip()

    # map long reads to contigs
    contigSam = contigFile.replace(".fasta", ".sam")
    command = ["minimap2", "-t", str(num_threads)]
    if read_set.platform == 'nanopore':
        command.extend(["-x", "map-ont"])
    elif read_set.platform == 'pacbio':
        command.extend(["-x", "map-pb"])
    command.extend(["-a", "-o", contigSam, contigFile, read_set.files[0]])
    tempTime = time()
    LOG.write(' '.join(command)+"\n")
    report['command'] = ' '.join(command)
    return_code = subprocess.call(command, shell=False, stderr=subprocess.DEVNULL)
    if return_code != 0:
        LOG.write("minimap2 map return code = %d, time = %d seconds\n"%(return_code, time() - tempTime))
        return None
    report['sam'] = contigSam

    if outformat == 'bam':
        (contigBam, samtools_version) = convertSamToBam(contigSam, num_threads, save_sam=save_intermediate)
        report['bam'] = contigBam
        report['samtools_version'] = samtools_version
    return report
            
def runBowtie(contigFile, read_set, num_threads, outformat='bam', save_intermediate=False):
    """
    index contigsFile, then run bowtie2, then convert sam file to pos-sorted bam and index
    """
    LOG.write("runBowtie(%s, %s, %s, %s) %s\n"%(contigFile, read_set.files[0], str(num_threads), outformat, strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    report = {}
    report['program'] = 'bowtie2'
    command = ["bowtie2", '--version']
    proc = subprocess.run(command, shell=False, text=True, capture_output=True)
    m = re.search("version\s+(\S+)", proc.stdout + proc.stderr)
    if m:
        report['version'] = m.group(1)

    if os.path.exists("bowtie_index_dir"):
        # delete dir and all files there
        shutil.rmtree('bowtie_index_dir')
    os.mkdir("bowtie_index_dir")
    command = ["bowtie2-build", "--threads", str(num_threads), contigFile, 'bowtie_index_dir/'+contigFile]
    LOG.write("executing: "+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    LOG.write("bowtie2-build return code = %d\n"%return_code)
    if return_code != 0:
        LOG.write("bowtie2-build failed\n")
        return None

    fastqBase = read_set.files[0]
    fastqBase = re.sub(r"\..*", "", fastqBase)
    samFile = contigFile+"_"+fastqBase+".sam"

    command = ["bowtie2", "-p", str(num_threads)]
    command.extend(["-x", 'bowtie_index_dir/'+contigFile])
    if len(read_set.files) > 1:
        #sys.stderr.write("we have a pair of read files\n")
        command.extend(('-1', read_set.files[0], '-2', read_set.files[1]))
    else:
        #sys.stderr.write("we have a list with a single read file\n")
        command.extend(('-U', read_set.files[0]))
    if read_set.format == 'fasta':
        command.append('-f')
    command.extend(('-S', samFile))
    LOG.write(" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False) #, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    LOG.write("bowtie2 return code = %d\n"%return_code)
    #shutil.rmtree('bowtie_index_dir')
    report['command'] = ' '.join(command)
    report['sam'] = samFile
    if not save_intermediate:
        shutil.rmtree('bowtie_index_dir')# delete dir and all files there
    if outformat == 'bam':
        (contigsBam, samtools_version) = convertSamToBam(samFile, num_threads,  save_sam=save_intermediate)
        LOG.write('got bam file %s\n'%contigsBam)
        report['bam'] = contigsBam
        report['samtools_version'] = samtools_version
    return report

def convertSamToBam(samFile, num_threads, save_sam=False):
    #convert format to bam and index
    LOG.write("convertSamToBam(%s, %s)\n"%(samFile, str(num_threads)))
    tempTime = time()
    command = ["samtools"]
    proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE)
    proc.wait()
    version_text = proc.stderr.read().decode()
    samtools_version = ''
    for line in version_text.splitlines():
        if 'Version' in line:
            samtools_version = line.strip()

    sortThreads = max(int(num_threads/2), 1)
    samFilePrefix = re.sub(".sam", "", samFile, re.IGNORECASE)
    command = ["samtools", "view", "-bS", "-@", str(sortThreads), "-o", samFilePrefix+"_unsorted.bam", samFile]
    LOG.write(" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    #os.remove(samFile) #save a little space
    if return_code != 0:
        comment = "samtools view returned %d"%return_code
        LOG.write(comment+"\n")
        return None
    if not save_sam:
        os.unlink(samFile)

    LOG.flush()

    bamFileSorted = samFilePrefix+".bam" 
    command = ["samtools", "sort", "-@", str(sortThreads), "-o", bamFileSorted, samFilePrefix+"_unsorted.bam"]
    if "Version: 0.1.19" in samtools_version:
        # different invocation for this older version
        command = ["samtools", "sort", "-@", str(sortThreads), samFilePrefix+"_unsorted.bam", samFilePrefix]
    return_code = subprocess.check_call(command, shell=False, stderr=LOG)

    if return_code != 0:
        comment = "samtools sort returned %d, convertSamToBam failed"%return_code
        LOG.write(comment+"\n")
        return None
    LOG.write("bamFileSorted = "+bamFileSorted+"\n")
    os.remove(samFilePrefix+"_unsorted.bam")
    if not os.path.exists(bamFileSorted):
        comment = "{0} not found, sorting bamfile failed, convertSamToBam failed\n".format(bamFileSorted)
        LOG.write(comment+"\n")
        return None
    if not os.path.getsize(bamFileSorted):
        comment = "{0} of size zero, sorting bamfile failed, convertSamToBam failed\n".format(bamFileSorted)
        LOG.write(comment+"\n")
        return None
    LOG.write("samtools sort return code=%d, time=%d, size of %s is %d\n"%(return_code, time()-tempTime, bamFileSorted, os.path.getsize(bamFileSorted)))

    command = ["samtools", "index", bamFileSorted]
    LOG.write("executing: "+" ".join(command)+"\n")
    return_code = subprocess.call(command, shell=False, stderr=LOG)
    #LOG.write("samtools index return code = %d\n"%return_code)
    return (bamFileSorted, samtools_version)

def runRacon(contigFile, read_set, num_threads, save_intermediate=False):
    """
    Polish (correct) sequence of assembled contigs by comparing to the original long-read sequences
    Run racon on reads, read-to-contig-sam, contigs. Generate polished contigs.
    Return name of polished contigs.
    """
    LOG.write('runRacon(%s, %s, %s)\n'%(contigFile, read_set.files[0], str(num_threads)))
    report = {"input_contigs":contigFile, "reads": read_set.files[0], "program": "racon"}
    command = ["racon", "--version"]
    proc = subprocess.Popen(command, shell=False, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    proc.wait()
    version_text = proc.stderr.read().decode()
    if not version_text:
        version_text = proc.stdout.read().decode()
    racon_version = version_text.strip()
    report['program'] = 'racon'
    report['version'] = racon_version

    minimap_report = runMinimap(contigFile, read_set, num_threads, outformat='sam', save_intermediate=save_intermediate)
    report['minimap report'] = minimap_report
    readsToContigsSam = minimap_report['sam']
    if not readsToContigsSam:
        comment = "runMinimap failed to generate sam file, exiting runRacon"
        LOG.write(comment + "\n")
        return None

    command = ["racon", "-t", str(num_threads), "-u"]
    averageQuality = read_set.avg_quality
    if averageQuality:
        command.extend(['-q', "{:.2f}".format(averageQuality * 0.5)])
    command.extend([ read_set.files[0], readsToContigsSam, contigFile])
    LOG.write("racon command: \n"+' '.join(command)+"\n")
    raconContigs = contigFile.replace(".fasta", ".racon.fasta")
    report['racon command'] = ' '.join(command)

    raconStartTime = time()
    with open(raconContigs, 'w') as raconOut:
        FNULL = open(os.devnull, 'w') # send stdout to dev/null
        return_code = subprocess.call(command, shell=False, stderr=FNULL, stdout=raconOut)
    LOG.write("racon return code = %d, time = %d seconds\n"%(return_code, time()-raconStartTime))
    report['seconds'] = time()-raconStartTime
    if return_code != 0:
        return None
    if not save_intermediate:
        os.remove(readsToContigsSam)
    raconContigSize = os.path.getsize(raconContigs)
    LOG.write("size of raconContigs: %d\n"%raconContigSize)
    if raconContigSize < 10:
        return None
    if re.search("racon.racon.fasta", raconContigs):
        shorterFileName = re.sub("racon.racon.fasta", "racon.fasta", raconContigs)
        report['original_name'] = raconContigs
        shutil.move(raconContigs, shorterFileName)
        LOG.write("renaming {} to {}\n".format(raconContigs, shorterFileName))
        raconContigs = shorterFileName
    report['output'] = raconContigs
    comment = "racon, input %s, output %s"%(contigFile, raconContigs)
    LOG.write(comment+"\n")
    #os.remove(contigFile) #delete old one
    return report

def runPilon(contigFile, read_set, num_threads, pilon_jar=None, save_intermediate=False):
    """ 
    polish contigs with short reads (illumina or iontorrent)
    first map reads to contigs with bowtie
    """
    LOG.write("runPilon starting Time = %s\n"%(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))))
    LOG.write("runPilon(%s, %s, %s)\n"%(contigFile, read_set.files[0], str(num_threads)))
    bowtieReport = runBowtie(contigFile, read_set, num_threads, outformat='bam', save_intermediate=save_intermediate)
    if not bowtieReport:
        sys.stderr.write("Problem: runBowtie failed")
        return None
    if not (os.path.exists(bowtieReport['bam'])):
        sys.stderr.write("Problem: runBowtie failed to return expected bamfile {}\n".format(bowtieReport['bam']))
        return None
    pilon_start_time = time()
    if (pilon_jar and os.path.exists(pilon_jar)):
        command = ['java', '-Xmx32G', '-jar', pilon_jar, '--genome', contigFile]
    else:
        command = ['pilon', '--genome', contigFile]
    if len(read_set.files) > 1:
        command.extend(('--frags', bowtieReport['bam']))
    else:
        command.extend(('--unpaired', bowtieReport['bam']))
    pilonPrefix = contigFile.replace(".fasta", "")
    m = re.match(".*pilon_(\d+)", pilonPrefix)
    if m :
        level = int(m.group(1))
        pilonPrefix = re.sub("pilon_"+m.group(1), "pilon_{}".format(level+1), pilonPrefix)
    else:
        pilonPrefix += "_pilon_1"
    pilonContigs = pilonPrefix  # + ".fasta" this gets added by pilon
    command.extend(('--outdir', '.', '--output', pilonContigs, '--changes'))
    command.extend(('--threads', str(num_threads)))
    LOG.write(" ".join(command)+"\n")
    LOG.write("bamfile = {}, size={}\n".format(bowtieReport['bam'], os.path.getsize(bowtieReport['bam'])))
    return_code = subprocess.call(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    LOG.write("pilon return code = %d\n"%return_code)
    pilon_time = time() - pilon_start_time
    LOG.write("pilon duration = %d\n"%(pilon_time))
    LOG.flush()
    if not save_intermediate:
        os.remove(bowtieReport['bam'])
        if os.path.exists(bowtieReport['bam']+".bai"):
            os.remove(bowtieReport['bam']+".bai")
    if return_code != 0:
        return None
    pilon_changes = 0
    with open(pilonContigs+".changes") as CHANGES:
        pilon_changes = len(CHANGES.read().splitlines())
    if not save_intermediate:
        os.remove(pilonContigs+".changes")
    if (pilon_jar and os.path.exists(pilon_jar)):
        command = ["java", "-jar", pilon_jar, "--version"]
    else:
        command = ['pilon', '--version']
    proc = subprocess.run(command, shell=False, capture_output=True, text=True)
    pilon_version = proc.stdout

    report = {"input_contigs":contigFile, 
            "reads": ":".join(read_set.files), 
            "program": "pilon", 
            "version": pilon_version, 
            "output": pilonContigs+".fasta", 
            "num_changes": pilon_changes, 
            "seconds" : pilon_time,
            "bowtie_version": bowtieReport['version']}
    if 'samtools_version' in bowtieReport:
        report['samtools_version'] = bowtieReport['samtools_version']

    comment = "pilon, input %s, output %s, num_changes = %d"%(contigFile, pilonContigs, pilon_changes)
    LOG.write(comment+"\n")
    if not save_intermediate:
        os.remove(contigFile) # clean up old one
    return report 

def calcReadDepth(bamfiles):
    """ Return dict of contig_ids to tuple of (coverage, normalized_coverage) """
    LOG.write("calcReadDepth(%s)\n"%" ".join(bamfiles))
    readDepth = {}
    command = ["samtools", "depth"]
    if type(bamfiles) is str:
        command.append(bamfiles)
    else:
        command.extend(bamfiles)
    LOG.write("command = "+" ".join(command)+"\n")
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    depthData = proc.communicate()[0].decode()
    if 0:
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

    #LOG.write("total length for depth data = %d\n"%totalLength)
    #LOG.write("total depth = %.1f\n"%totalDepthSum)
    LOG.write("len(readDepth) = %d\n"%len(readDepth))
    totalMeanDepth = 0
    if totalLength > 0:
        totalMeanDepth = totalDepthSum/totalLength

    # calculate mean depth of contigs within "normal" boundary around overall mean
    lowerBound = totalMeanDepth * 0.5
    upperBound = totalMeanDepth * 2
    oneXSum = 0.0
    oneXLen = 0
    for c in readDepth:
        if c == 'average_coverage':
            continue
        meanDepth = readDepth[c][0]
        if meanDepth >= lowerBound and meanDepth <= upperBound:
            oneXSum += meanDepth * contigLength[c]
            oneXLen += contigLength[c]
    oneXDepth = 1
    if oneXLen > 0 and oneXSum > 0:
        oneXDepth = oneXSum/oneXLen # length-weighted average

    for c in readDepth:
        if c == 'average_coverage':
            continue
        meanDepth = readDepth[c][0]
        normalizedDepth = meanDepth / oneXDepth
        readDepth[c][1] = normalizedDepth
    return (totalMeanDepth, readDepth)

def write_html_report(htmlFile, read_list, assembly_data):
    LOG.write("writing html report to %s\n"%htmlFile)
    HTML = open(htmlFile, 'w')
    HTML.write("<!DOCTYPE html><html><head>\n")
    HTML.write("<link href=\"https://fonts.googleapis.com/css?family=Work+Sans:300,400,500,600,700,800,900\" rel=\"stylesheet\">\n")
    HTML.write("""
<style>
 body {
    font-family: 'Work Sans', sans-serif;
    color: #444;
    }
 header { padding-bottom: .5in; }
 section { margin-bottom: .5in; }

 a {
    text-decoration: none;
    color: #0d78ef;
    }
 a:hover { text-decoration: underline; }
 h2 { font-size: 1.45em; }
 h2, h3 { font-weight: 500; }
 h2 small { color: #888; font-size: 0.65em; }
 .pull-left { float: left; }
 .pull-right { float: right; }
 sup { display: inline-block; }
 .code {font-family: monospace; }

 /* tables */
 table {
     margin: .1in 0 .5in 0;
     border-collapse: collapse;
     page-break-inside: avoid;
     }
 .table-header {
     color: #fff;
     background: #196E9C;
     padding: 5px;
     text-align: left;
     }
 .table-header th { font-weight: 400; }
 .row-header { border-bottom: 1px solid #196E9C; font-weight: 600; }
 th, td { padding: 5px; }
 th { border-bottom: 2px solid #196E9C; }
 tr:last-child { border-bottom: 1px solid #196E9C; }

 .kv-table
     text-align: left;
     }

 table td.align-right {
     text-align: right;
     }

 .lg-table { width: 80%; }
 .med-table { width: 60%; }
 .sm-table { width: 40%; }
 .xs-table {width: 20%; }

 .center { margin-left: auto; margin-right: auto; }

 .title {
 padding: 17px 0 0 0px;
 font-size: 1.3em;
 color: #777;
 }

 """)
    
    program_version = {}
    HTML.write("</style></head><body><header>\n")
    HTML.write('<div class="report-info pull-right">\n')
    HTML.write('<span>Report Date:</span> '+strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+'<br>\n')
    HTML.write('</div></header>\n')
    HTML.write("<h1>Genome Assembly Report</h1>\n")

    if "bandage" in assembly_data and "plot" in assembly_data["bandage"]:
        bandagePlotFile = assembly_data["bandage"]["plot"]
        if os.path.exists(bandagePlotFile):
            HTML.write('<section>\n<h2>Bandage Plot</h2>\n<div class="a">\n')
            svg_text = open(bandagePlotFile).read()
            svg_text = re.sub(r'<svg width="[\d\.]+mm" height="[\d\.]+mm"', '<svg width="200mm" height="150mm"', svg_text)
            HTML.write(svg_text+"\n\n</div>\n")
            program_version[assembly_data["bandage"]['program']] = assembly_data["bandage"]["version"]
        HTML.write("</section>\n")
    elif 'bandage_missing_message' in assembly_data:
        HTML.write("<section>"+assembly_data['bandage_missing_message']+"</section>")

    HTML.write('<section>\n<h2>Assembly</h2>\n')
    HTML.write("""
    <table class="med-table kv-table">
        <thead class="table-header">
        <tr> <th colspan="2"> Assembly Process </th></tr></thead>
        <tbody>
        """)
    for key in ['program', 'version', 'duration', 'contigs file size']:
        HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, str(assembly_data[key])))
    HTML.write("<tr><td class='code'>command:</td><td>%s</td></tr>\n"%(str(assembly_data['command'])))
    HTML.write("</tbody></table>\n")
    HTML.write("</section>\n")
    program_version[assembly_data['program']] = assembly_data["version"]

    if 'polishing' in assembly_data and len(assembly_data['polishing']):
        HTML.write('<section>\n<h2>Polishing</h2>\n')
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> Polishing Rounds </th></tr></thead>
            <tbody>
            """)
        for iteration, info in enumerate(assembly_data['polishing']):
            if info:
                HTML.write("<tr class='row-header'><td>%s:</td><td>%d</td></tr>\n"%("Round", iteration+1))
                for key in ['program', 'reads', 'seconds']:
                    value = str(info[key])
                    if key == 'seconds':
                        value = "{:.2f}".format(info[key])
                    HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(key, value))
                HTML.write("<tr></tr>\n") # blank row
                program_version[info['program']] = info["version"]
                if 'bowtie2_version' in info:
                    program_version['bowtie2'] = info['bowtie2_version']
                if 'samtools_version' in info:
                    program_version['samtools'] = info['samtools_version']
        HTML.write("</tbody></table>\n")
        HTML.write("</section>\n")

    if "quast" in assembly_data:
        HTML.write("<section><h2>Quast Report</h2>\n")
        if 'report_text' in assembly_data["quast"]:
            HTML.write("<pre>\n")
            HTML.write(assembly_data["quast"]["report_text"])
            HTML.write("\n</pre>\n")
        HTML.write("</section>\n")
        program_version['quast'] = assembly_data['quast']["version"]
    
    if "contig_filtering" in assembly_data and assembly_data['contig_filtering']:
        HTML.write('<section>\n<h2>Filtering Contigs on Length and Coverage</h2>\n')
        HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th colspan="2"> Contig Filtering </th></tr></thead>
            <tbody>
            """)
        for key in sorted(assembly_data['contig_filtering']):
            HTML.write("<tr><td>{}:</td><td>{}</td></tr>\n".format(key, assembly_data['contig_filtering'][key])) 
        HTML.write("</tbody></table>\n")
        HTML.write("</section>\n")
        if 'bowtie2_version' in assembly_data['contig_filtering']:
            program_version['bowtie2'] = assembly_data['contig_filtering']['bowtie2_version']
        if 'samtools_version' in assembly_data['contig_filtering']:
            program_version['samtools'] = assembly_data['contig_filtering']['samtools_version']

    HTML.write('<div>\n<h2>Preprocessing of Reads</h2>\n')
    for read_set in read_list:
        read_set.writeHtml(HTML)
    HTML.write("</div>\n")

    HTML.write("<section><h2>Tools Used:</h2>\n")
    HTML.write("""
        <table class="med-table kv-table">
            <thead class="table-header">
            <tr> <th > Tool </th><th> Version</th></tr></thead>
            <tbody>
            """)
    for program in ReadLibrary.program_version:
        program_version[program] = ReadLibrary.program_version[program]
    for program in sorted(program_version):
        HTML.write("<tr><td>%s:</td><td>%s</td></tr>\n"%(program, str(program_version[program])))
    HTML.write("</table>\n")
    HTML.write("</section>\n")
    HTML.write("</html>\n")

    HTML.close()


def main():
    main_return_code = 1 # set to zero when we have an assembly
    START_TIME = time()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outputDirectory', '-d', default='p3_assembly_work')
    illumina_or_iontorrent = parser.add_mutually_exclusive_group()
    illumina_or_iontorrent.add_argument('--illumina', metavar='files', nargs='*', help='Illumina fastq[.gz] files or pairs; use ":" between end-pairs or  percent-sign between mate-pairs', required=False, default=[])
    illumina_or_iontorrent.add_argument('--iontorrent', metavar='files', nargs='*', help='list of IonTorrent[.gz] files or pairs, use : between paired-end-files', required=False, default=[])
    parser.add_argument('--pacbio', metavar='files', nargs='*', help='list of Pacific Biosciences fastq[.gz] files', required=False, default=[])
    parser.add_argument('--nanopore', metavar='files', nargs='*', help='list of Oxford Nanotech fastq[.gz] files', required=False, default=[])
    parser.add_argument('--fasta', metavar='files', nargs='*', help='list of fasta[.gz] files', required=False, default=[])
    parser.add_argument('--sra', metavar='files', nargs='*', help='list of SRA run accessions (e.g. SRR5070677), will be downloaded from NCBI', required=False)
    parser.add_argument('--anonymous_reads', metavar='files', nargs='*', help='unspecified read files, types automatically inferred.')
    parser.add_argument('--max_bases', type=int, default=DEFAULT_MAX_BASES, help='downsample reads if more than this (per read set).')
    parser.add_argument('--interleaved', nargs='*', help='list of fastq files which are interleaved pairs')
    parser.add_argument('--recipe', choices=['unicycler', 'flye', 'meta-flye', 'canu', 'spades', 'meta-spades', 'plasmid-spades', 'single-cell', 'rna-spades', 'megahit', 'auto', 'none'], help='assembler to use', default='auto')
    parser.add_argument('--contigs', metavar='fasta', help='perform polishing on existing assembly')
    #parser.add_argument('--only-assembler', action='store true', help='omit spades read error correction')
    
    parser.add_argument('--racon_iterations', type=int, default=0, help='number of times to run racon per long-read file', required=False)
    parser.add_argument('--pilon_iterations', type=int, default=0, help='number of times to run pilon per short-read file', required=False)
    parser.add_argument('--pilon_hours', type=float, default=6.0, help='maximum hours to run pilon', required=False)
    parser.add_argument('--prefix', default='', help='prefix for output files', required=False)
    parser.add_argument('--genome_size', metavar='base pairs', type=int, default=DEFAULT_GENOME_SIZE, help='genome size for canu or flye, eg 5000000', required=False)
    parser.add_argument('--min_contig_length', type=int, default=300, help='save contigs of this length or longer', required=False)
    parser.add_argument('--min_contig_coverage', type=float, default=5, help='save contigs of this coverage or deeper', required=False)
    #parser.add_argument('--fasta', nargs='*', help='list of fasta files "," between libraries', required=False)
    parser.add_argument('--trusted_contigs', help='for SPAdes, same-species contigs known to be good', required=False)
    parser.add_argument('--no_pilon', action='store_true', help='for unicycler', required=False)
    parser.add_argument('--untrusted_contigs', help='for SPAdes, same-species contigs used gap closure and repeat resolution', required=False)
    parser.add_argument('-t', '--threads', metavar='cores', type=int, default=8)
    parser.add_argument('-m', '--memory', metavar='GB', type=int, default=DEFAULT_MEMORY_GIGABYTES, help='RAM limit in Gb')
    parser.add_argument('--trim', action='store_true', help='trim reads with trim_galore at default settings')
    parser.add_argument('--normalize', action='store_true', help='normalize read depth to target depth with BBNorm')
    parser.add_argument('--filtlong', action='store_true', help='filter long reads to target depth using filtlong')
    parser.add_argument('--target_depth', type=int, default=200, help='downsample to this approx. coverage')
    parser.add_argument('--pilon_jar', help='path to pilon executable or jar', default=None)
    parser.add_argument('--canu_exec', default="canu", help='path to canu executable (def "canu")')
    parser.add_argument('--spades_for_unicycler', help='path to spades.py suitable for unicycler')
    parser.add_argument('--bandage', action='store_true', help='generate image of assembly path using Bandage')
    parser.add_argument('--maxContigsForBandage', type=int, default=100, help='only generage Bandage for this or fewer contigs')
    parser.add_argument('--params_json', help='JSON file with additional information.')
    parser.add_argument('--path_prefix', '--path-prefix', help="Add the given directories to the PATH", nargs='*', required=False)
    parser.add_argument('-v', '--verbose', action='store_true', help="Output more progress information.", required=False)
    parser.add_argument('--debug', action='store_true', help="Save intermediate files.", required=False)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()

    if args.prefix:
        args.prefix = re.sub("[^\w\d\-_.]", "_", args.prefix)
        if not args.prefix.endswith("_"):
            args.prefix += "_"
    WORK_DIR = args.outputDirectory
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
    #logfileName = os.path.join(DETAILS_DIR, args.prefix + "p3_assembly.log")
 
    if args.path_prefix:
        os.environ["PATH"] = ":".join(args.path_prefix) + ":" + os.environ["PATH"]

    global LOG 
    #sys.stderr.write("logging to "+logfileName+"\n")
    #LOG = open(logfileName, 'w') 
    LOG = sys.stderr
    LOG.write("starting %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(START_TIME))+"\n")
    LOG.write("args= "+str(args)+"\n\n")
    LOG.write("Work directory is "+WORK_DIR+"\n\n")
    LOG.write("Final output will be saved to "+SAVE_DIR+"\n\n")
    LOG.write("Detailed output will be saved to "+DETAILS_DIR+"\n\n")

    ReadLibrary.NUM_THREADS = args.threads
    ReadLibrary.MEMORY = args.memory  # in GB
    ReadLibrary.LOG = sys.stderr

    read_list = []
    if args.anonymous_reads:
        for item in args.anonymous_reads:
            if ':' in item:
                read_pair = item.split(':')
                readLib = ReadLibrary(read_pair, work_dir=WORK_DIR)
            else:
                readLib = ReadLibrary(item, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.illumina:
        platform = 'illumina'
        for item in args.illumina:
            if ':' in item:
                read_pair = item.split(':')
                readLib = ReadLibrary(read_pair, platform=platform, work_dir=WORK_DIR)
            else:
                interleaved = args.interleaved and item in args.interleaved
                readLib = ReadLibrary(item, platform=platform, interleaved=interleaved, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.iontorrent:
        platform = 'iontorrent'
        for item in args.iontorrent:
            interleaved = args.interleaved and item in args.interleaved
            readLib = ReadLibrary(item, platform=platform, interleaved=interleaved, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.pacbio:
        platform = 'pacbio'
        for item in args.pacbio:
            readLib = ReadLibrary(item, platform=platform, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.nanopore:
        platform = 'nanopore'
        for item in args.nanopore:
            readLib = ReadLibrary(item, platform=platform, work_dir=WORK_DIR)
            read_list.append(readLib)

    if args.fasta:
        for item in args.fasta:
            if ':' in item:
                read_pair = item.split(':')
                readLib = ReadLibrary(read_pair, work_dir=WORK_DIR)
            else:
                readLib = ReadLibrary(item, work_dir=WORK_DIR)
            read_list.append(readLib)

    # move into working directory so that all files are local
    os.chdir(WORK_DIR)

    for read_set in read_list:
        read_set.study_reads()

    constrain_total_bases(read_list, args.max_bases)

    any_short_fasta = False
    short_reads = []
    long_reads = []
    for read_set in read_list:
        if read_set.length_class == 'short':
            short_reads.append(read_set)
            if read_set.format == 'fasta':
                any_short_fasta = True
        else:
            long_reads.append(read_set)

    if args.trim:
        for read_set in short_reads:
            if read_set.format == 'fastq': # TrimGalore only works on short fastq reads
                read_set.trim_short_reads()

    if args.normalize:
        LOG.write(f"args.normalize is set\n")
        for read_set in short_reads:
            read_set.normalize_read_depth(target_depth = int(args.target_depth/len(read_list)))

    if args.filtlong and len(long_reads) and args.target_depth:
        LOG.write(f"args.filtlong is set\n")
        illumina_reference = None
        for read_set in short_reads:
            if len(read_set.files) > 1:
                illumina_reference = read_set
                LOG.write(f"selected illumina reference for filtlong: {illumina_reference}\n")
                break
        for read_set in long_reads:
            target_bases = int(args.genome_size * args.target_depth / len(read_list))
            if target_bases < read_set.num_bases:
                read_set.filter_long_reads(target_bases, illumina_reference)

    if args.recipe == "auto":
        #now must decide which assembler to use
        LOG.write("translate auto into specific recipe\n")
        LOG.write("number of short reads = {}\n".format(len(short_reads)))
        LOG.write("number of long reads = {}\n".format(len(long_reads)))
        if len(short_reads):
            if any_short_fasta:
                args.recipe = "spades"
                LOG.write("auto recipe selecting spades due to presence of short fasta read data.\n")
            else:
                args.recipe = "unicycler"
                LOG.write("auto recipe selecting unicycler due to presence of short fastq read data.\n")
        elif len(long_reads):
            args.recipe = "flye"
            LOG.write("auto recipe selecting flye due to presence of long and absence of short read data.\n")
        else:
            comment = "auto recipe failed to find short or long reads\n"
            LOG.write(comment)
            raise Exception(comment)

    if ("meta" in args.recipe) or (args.recipe == 'megahit'): # meta-spades or meta-flye or megahit
        args.pilon_iterations = 0
        args.racon_iterations = 0
        comment = "Because recipe is for metagenome, turning pilon and racon iterations off."
        LOG.write(comment+"\n")

    contigs = ""
    startTime = time()
    assembly_data = None
    if args.recipe == "unicycler":
        spades_exec = None
        if (args.spades_for_unicycler):
            spades_exec = args.spades_for_unicycler
        assembly_data = runUnicycler(read_list, threads=args.threads, min_contig_length=args.min_contig_length, prefix=args.prefix, spades_exec=spades_exec )

    elif args.recipe == "canu":
        assembly_data = runCanu(read_list, canu_exec=args.canu_exec, threads=args.threads, genome_size=args.genome_size, memory=args.memory, prefix=args.prefix)

    elif "flye" in args.recipe:
        assembly_data = runFlye(read_list, recipe=args.recipe, threads=args.threads, genome_size=args.genome_size, prefix=args.prefix)

    elif "spades" in args.recipe or args.recipe == "single-cell":
        assembly_data = runSpades(read_list, prefix=args.prefix, recipe=args.recipe, threads=args.threads, memoryGB=args.memory)
    elif args.recipe == 'megahit':
        assembly_data = runMegahit(read_list, prefix=args.prefix, threads=args.threads, memoryGB=args.memory)
    elif args.recipe == 'none':
        LOG.write("recipe specified as 'none', no assembly will be performed\n")
    else:
        LOG.write("cannot interpret args.recipe: "+args.recipe)

    elapsedTime = time() - startTime
    elapsedHumanReadable = ""
    if elapsedTime < 60:
        elapsedHumanReadable = "%.1f minutes"%(elapsedTime/60.0)
    elif elapsedTime < 3600:
        elapsedHumanReadable = "%.2f hours"%(elapsedTime/3600.0)
    else:
        elapsedHumanReadable = "%.1f hours"%(elapsedTime/3600.0)
    elapsedHumanReadable += f" ({int(elapsedTime)} seconds)"

    if os.path.exists('contigs.fasta'):
        contigs = 'contigs.fasta'
        assembly_data["contigs file size"] = os.path.getsize("contigs.fasta")
        assembly_data['duration'] = elapsedHumanReadable
        LOG.write("size of contigs file is %d\n"%os.path.getsize(contigs))
    else:
        comment = "assembly failed to generate contigs"
        LOG.write(comment+"\n")
        sys.exit(1) # signal job failure

    polishing_records = []
    if args.racon_iterations and contigs:
        # now run racon with each long-read file
            for longReadSet in long_reads:
                for i in range(0, args.racon_iterations):
                    LOG.write("runRacon on {}, {}, platform={}, round={}\n".format(contigs, longReadSet.files[0], longReadSet.platform, i))
                    raconReport = runRacon(contigs, longReadSet, args.threads, save_intermediate=args.debug)
                    raconContigFile = ''
                    if raconReport:
                        if 'output' in raconReport:
                            raconContigFile = raconReport['output']
                        polishing_records.append(raconReport)
                        LOG.write("racon output = {}, report={}\n".format(raconContigFile, raconReport))
                        if os.path.exists(raconContigFile):
                            contigs = raconContigFile
                            sys.stderr.write("contigs file is now {}\n".format(contigs))
                        else:
                            break # break out of iterating racon_iterations, go to next long-read file if any
                    else:
                        LOG.write("runRacon failed to return expected data")
                        break
        
    if args.pilon_iterations and contigs:
        pilon_end_time = time() + args.pilon_hours * 60 * 60
        # now run pilon with each short-read file
        for read_set in short_reads:
            fastqFile = read_set.files[0] # use only read_1 if paired
            for iteration in range(0, args.pilon_iterations):
                if time() > pilon_end_time:
                    LOG.write("Time expended on pilon exceeds allocation of {} hours. Omitting further pilon runs.".format(args.pilon_hours))
                    break

                LOG.write("runPilon(%s, %s, ...)\n"%(contigs, fastqFile))
                pilonReport = runPilon(contigs, read_set, args.threads, args.pilon_jar, save_intermediate=args.debug)
                if pilonReport:
                    polishing_records.append(pilonReport)
                    pilonContigFile = pilonReport['output']
                    if pilonContigFile is not None and os.path.exists(pilonContigFile):
                        contigs = pilonContigFile
                        sys.stderr.write("contigs file is now {}\n".format(contigs))
                    else:
                        sys.stderr.write("expected contigs file {} does not exist\n".format(contigs))
                        #break
                    # check number of changes in most recent round, quit if zero
                    if ('num_changes' in pilonReport) and (pilonReport['num_changes'] == 0):
                        break
                else:
                    LOG.write("runPilon failed to return expected data")
                    break
    assembly_data['polishing'] = polishing_records
            
    if os.path.exists(contigs) and os.path.getsize(contigs):
        filterReport = filterContigsByLengthAndCoverage(contigs, read_list, args.min_contig_length, args.min_contig_coverage, args.threads, args.prefix)
        assembly_data['contig_filtering'] = filterReport
        if 'good contigs file' in filterReport:
            contigs = filterReport['good contigs file']
        if contigs and os.path.getsize(contigs):
            quastReport = runQuast(contigs, args.min_contig_length, args.threads)
            if quastReport:
                assembly_data['quast'] = quastReport
                if 'file_path' in quastReport:
                    shutil.move(quastReport["file_path"], os.path.join(DETAILS_DIR, args.prefix+"quast_report.txt"))
            shutil.move(contigs, os.path.join(SAVE_DIR, args.prefix+"contigs.fasta"))

        gfaFile = os.path.join(DETAILS_DIR, args.prefix+"assembly_graph.gfa")
        if os.path.exists(gfaFile) and os.path.getsize(gfaFile):
            if 'total contigs' in filterReport:
                if int(filterReport['total contigs']) <= args.maxContigsForBandage:
                    bandageReport = runBandage(gfaFile)
                    assembly_data['bandage'] = bandageReport
                else:
                    assembly_data['bandage_missing_message'] = f"Bandage plot not generated for assemblies with over {args.maxContigsForBandage} contigs."


    else:
        LOG.write("contigs not generated")
        sys.exit(1)

    with open(os.path.join(DETAILS_DIR, args.prefix+"run_details.json"), "w") as fp:
        try:
            json.dump(assembly_data, fp, indent=2, sort_keys=True)
        except UnicodeDecodeError as ude:
            LOG.write("Problem writing assembly_data to json: "+str(ude)+"\n")


    htmlFile = os.path.join(SAVE_DIR, args.prefix+"AssemblyReport.html")
    write_html_report(htmlFile, read_list, assembly_data)
    LOG.write("done with %s\n"%sys.argv[0])
    LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
    LOG.write("Total time in hours = %d\t"%((time() - START_TIME)/3600))
    LOG.write("Total time in seconds = %d\n"%((time() - START_TIME)))
    LOG.close()
    return main_return_code

if __name__ == "__main__":
    main()
