import csv
import os
import json
import yaml
import numpy as np
from snakemake.logging import logger
import re

#shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")

configfile:config["all"]
FILES = json.load(open(config['SAMPLES_JSON']))
WORKDIR = os.path.abspath(config["OUTPUT_DIR"])
PROJECT_NAME = config['PROJECT_NAME']


###########################################################################
#################### Defining samples, cases, controls ####################
###########################################################################

SAMPLES_NAMES = sorted(FILES.keys())


def getChIP(wildcards):
	return FILES[wildcards.sample]['ChIP']

def getInput(wildcards):
	return FILES[wildcards.sample]['Input']

def getChIP_tagSize(wildcards):
	return FILES[wildcards.sample]['ChIP_tagSize']

def getInput_tagSize(wildcards):
	return FILES[wildcards.sample]['Input_tagSize']

sampleList = FILES.keys()
#wildcard_constraints:
#    sample = sampleList

rule all:
    input:
        expand("peaks/{sample}_peaks_q5_peaks.broakPeak",sample=sampleList)
		
rule filterDupInput:
	input:getInput
	output:bed = "tmpBdg/{sample}_inputfilterDup.bed",
		metric = "tmpBdg/{sample}_inputfilterDup.metrics.txt"
	shell:
		"""
        macs2 filterdup -i {input} --keep-dup 1 -o {output.bed} &> {output.metric}
        """

rule catchInput:
	input:"tmpBdg/{sample}_inputfilterDup.metrics.txt"
	output:"tmpBdg/{sample}_inputfilterDup.fragmentNb.txt"
	shell:
		"""
        grep "tags after filtering in alignment file: " {input} | awk -F " " '{{print $14}}' > {output}
        """

rule collectInsertSizeMetrics:
	input: getChIP
	output:"tmpInsertSizeMet/{sample}_insertSizeMetrics.txt",
	shell:
		"""
		set +o pipefail
		samtools view -h {input} | perl ~/Documents/tools/ChIPSeq/sam_insert-size.pl -i - -p -d -f -min 50 -max 500 -n 2000000 -xlim_i 350 -xlim_r 200 &> {output}
		set -o pipefail
		"""

rule catchMedFragLength:
	input:"tmpInsertSizeMet/{sample}_insertSizeMetrics.txt"
	output:"tmpInsertSizeMet/{sample}_medFragLength.txt"
	shell:
		"""
        grep "Insert size median: " {input} | awk -F ": " '{{printf "%d\\n%d",$2/2,$2}}' > {output}
        """

rule filterDupChIP:
	input: getChIP
	output: bed = temp("tmpBdg/{sample}_filterDup.bed"),
		metric = temp("tmpBdg/{sample}_filterDup.metrics.txt")
	shell:
		"""
		macs2 filterdup -i {input} --keep-dup 1 -o {output.bed} -f BAMPE &> {output.metric}
		"""

rule catchChIP:
	input:"tmpBdg/{sample}_filterDup.metrics.txt"
	output:"tmpBdg/{sample}_filterDup.fragmentNb.txt"
	shell:
		"""
        grep "tags after filtering in alignment file: " {input} | awk -F" " '{{print $14}}' > {output}
        """

rule predictd:
	input: "tmpBdg/{sample}_filterDup.bed"
	output: temp("tmpBdg/{sample}_predictd.metrics.txt")
	shell:
		"""macs2 predictd -i {input} -g hs -m 5 50 &> {output}"""

# ChIP pileup
rule pileup_ChIP:
	input:"tmpBdg/{sample}_filterDup.bed"
	output:temp("tmpBdg/{sample}_filterDup.pileup.bdg")
	shell:
		"""
        macs2 pileup -i {input} -o {output} -f BEDPE
        """

# d background
rule pileup_fregLeng_bg:
	input:
		bed="tmpBdg/{sample}_inputfilterDup.bed",
		FragFile="tmpInsertSizeMet/{sample}_medFragLength.txt"
	output: temp("tmpBdg/{sample}_d.bg.bdg")
	run:
		f=open("tmpInsertSizeMet/"+wildcards.sample+"_medFragLength.txt")
		half=f.readline().strip()
		f.close()
		shell("macs2 pileup -i {input.bed} -f BED -B --extsize {half} -o {output}")

# The slocal background building
rule pileup_1k_bg:
	input:
		bed="tmpBdg/{sample}_inputfilterDup.bed"
	output: temp("tmpBdg/{sample}_1k.bg.bdg")
	shell:
		"""macs2 pileup -f BED -i {input.bed} -B --extsize 500 -o {output}"""
rule opt_1k_bg:
	input:
		bdg="tmpBdg/{sample}_1k.bg.bdg",
		frag="tmpInsertSizeMet/{sample}_medFragLength.txt"
	output: temp("tmpBdg/{sample}_1k.bg.norm.bdg")
	run:
		f=open("tmpInsertSizeMet/"+wildcards.sample+"_medFragLength.txt")
		frag=int(f.readlines()[1].strip())/1000
		f.close()
		shell("macs2 bdgopt -i {input.bdg} -m multiply -p {frag} -o {output}")

# The llocal background building
rule pileup_10k_bg:
	input:
		bed="tmpBdg/{sample}_inputfilterDup.bed"
	output: temp("tmpBdg/{sample}_10k.bg.bdg")
	shell:
		"""macs2 pileup -f BED -i {input.bed} -B --extsize 5000 -o {output}"""
rule opt_10k_bg:
	input:
		bdg="tmpBdg/{sample}_10k.bg.bdg",
		frag="tmpInsertSizeMet/{sample}_medFragLength.txt"
	output: temp("tmpBdg/{sample}_10k.bg.norm.bdg")
	run:
		f=open("tmpInsertSizeMet/"+wildcards.sample+"_medFragLength.txt")
		frag=int(f.readlines()[1].strip())/10000
		f.close()
		shell("macs2 bdgopt -i {input.bdg} -m multiply -p {frag} -o {output}")


# MAX signal between slocal and llocal
rule max_Slocal_Llocal:
	input: slocal="tmpBdg/{sample}_1k.bg.norm.bdg",
		llocal="tmpBdg/{sample}_10k.bg.norm.bdg"
	output: temp("tmpBdg/{sample}_1k.10k.bg.norm.bdg")
	shell:
		"""macs2 bdgcmp -m max -t {input.slocal} -c {input.llocal} -o {output}"""

# MAX signal between slocal_llocal and d background
rule max_local_dFrag:
	input: local="tmpBdg/{sample}_1k.10k.bg.norm.bdg",
		d_bdg="tmpBdg/{sample}_d.bg.bdg"
	output: temp("tmpBdg/{sample}_d.1k.10k.bg.norm.bdg")
	shell:
		"""macs2 bdgcmp -m max -t {input.local} -c {input.d_bdg} -o {output}"""

# MAX signal between local_d and genome background
rule max_local_genome:
	input: local="tmpBdg/{sample}_d.1k.10k.bg.norm.bdg",
		genomeBdg="tmpInsertSizeMet/{sample}_medFragLength.txt",
		fragLen="tmpBdg/{sample}_inputfilterDup.fragmentNb.txt"
	output: temp("tmpBdg/{sample}_local.bias.raw.bdg")
	run:
		f=open("tmpInsertSizeMet/"+wildcards.sample+"_medFragLength.txt")
		frag=int(f.readlines()[1].strip())
		f.close()
		f=open("tmpBdg/"+wildcards.sample+"_inputfilterDup.fragmentNb.txt")
		ctlReadNb=int(f.readline().strip())
		f.close()
		genomeScale=float((frag*ctlReadNb)/2700000000)
		shell("macs2 bdgopt -m max -i {input.local} -p {genomeScale} -o {output}")

# Scale ChIP to control
rule scaleByLibSize:
	input:
		rawBias="tmpBdg/{sample}_local.bias.raw.bdg",
		Ctl="tmpBdg/{sample}_inputfilterDup.fragmentNb.txt",
		chip="tmpBdg/{sample}_filterDup.fragmentNb.txt"
	output:"tmpBdg/{sample}_local.lambda.bdg"
	run:
		f=open("tmpBdg/"+wildcards.sample+"_inputfilterDup.fragmentNb.txt")
		inputTot=int(f.readline().strip())
		f.close()
		f=open("tmpBdg/"+wildcards.sample+"_filterDup.fragmentNb.txt")
		chipTot=int(f.readline().strip())
		f.close()
		scale=float(inputTot/chipTot)
		shell("macs2 bdgopt -i {input.rawBias} -m multiply -p {scale} -o {output}")

# Compare ChIP and Control
rule compareChipControl:
	input:
		chip="tmpBdg/{sample}_filterDup.pileup.bdg",
		bias="tmpBdg/{sample}_local.lambda.bdg"
	output: "tmpBdg/{sample}_qvalue.bdg"
	shell:
		"""macs2 bdgcmp -t {input.chip} -c {input.bias} -m qpois -o {output}"""

rule callpeak:
	input:
		fragMet="tmpInsertSizeMet/{sample}_medFragLength.txt",
		chip="tmpBdg/{sample}_qvalue.bdg"
	output:"peaks/{sample}_peaks_q5_peaks.broakPeak"
	run:
		f=open("tmpInsertSizeMet/"+wildcards.sample+"_medFragLength.txt")
		frag=f.readlines()[1].strip()
		f.close()
		shell("macs2 bdgbroadcall -i {input.chip} -c 5 -l {frag} -g 101 -o {output}")
