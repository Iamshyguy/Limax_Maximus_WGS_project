import argparse
import subprocess
import os
import shutil

#set up parser
parser = argparse.ArgumentParser(
                    prog='assembly',
                    description='Script that performs a whole genome assembly using long read (nanopore or pacbio), short read (illumina) and hi-c (illumina) data. The script first assembles the long read data, polishes it using the short read data and then sorts it into chromosomes using the hi-c data',
                    epilog='')

parser.add_argument('long_input',help="path to the long read data fastq")
parser.add_argument('short_input1',help="path to the forward short read data fastq (please ensure clean_lowqual from clean_illumina is installed in your path)")
parser.add_argument('short_input2',help="path to the reverse short read data fastq (please ensure clean_lowqual from clean_illumina is installed in your path)")
parser.add_argument('hic_input1',help="path to the forward hi-c read data fastq.gz (please ensure these are gzipped to avoid juicer errors)")
parser.add_argument('hic_input2',help="path to the reverse hi-c read data fastq.gz (please ensure these are gzipped to avoid juicer errors)")
parser.add_argument('outpath',help="path to the output directory (will create one if the directory is not yet present, if present, will ask if user wants to remove or rename the directory)")
parser.add_argument('-p', '--pacbio', action='store_true', help="utilize pacbio hifi data instead of nanopore")
parser.add_argument('-s', '--snail', help="the name of the snail genome being assembled, default = lmaximus" )
parser.add_argument('-t', '--threads', help="the amount of threads you'd like to use (default = 94)")
parser.add_argument('-j', '--juicerdir', help="location of the juicer CPU directory (default = /media/data2/juicer-1.6/CPU)")
parser.add_argument('-d', '--dna3ddir', help="location of the 3d-dna directory (default = /media/data2/3d-dna-master_altered_fixed)")
parser.add_argument('-pj', '--pilonjar', help="location of the pilon jar (default = /home/milan/pilon-1.24.jar )")
parser.add_argument('-sf', '--start_from', help="the step of the assembly to start from, to see the list of steps use -sl")
parser.add_argument('-sl', '--stepslist', action='store_true',help="Showcase a list of steps for the -sf flag")
parser.add_argument('-ha', '--hifiasm_args', help="additional hifiasm arguments")
parser.add_argument('-pa', '--pilon_args', help="additional pilon arguments")
parser.add_argument('-ja', '--juicer_args', help="additional hifiasm arguments")
parser.add_argument('-da', '--dna3d_args', help="additional hifiasm arguments")
parser.add_argument('-sj', '--stop_at_juicer', action='store_true', help="Stop the pipeline just before the juicer step. Due to the way juicer is coded, only 1 run of the pipeline can be performed at a time unless this flag is used and the juicer and 3DDNA steps are then ran seperately")
parser.add_argument('-ns', '--no_short_reads', action='store_true', help="Skip steps involving illumina genomic short reads (filtering, mapping and pilon correction)")



args = parser.parse_args()

ns=args.no_short_reads
sj=args.stop_at_juicer
longr=args.long_input
shortr1=args.short_input1
shortr2=args.short_input2
hicr1=args.hic_input1
hicr2=args.hic_input2
pac=args.pacbio
snail=args.snail.replace(" ","")
outfile=args.outpath
threads=args.threads
juicerdir=args.juicerdir
dnadir=args.dna3ddir
pilonjar=args.pilonjar
hifiasm_args=args.hifiasm_args
pilon_args=args.pilon_args
juicer_args=args.juicer_args
dna3d_args=args.dna3d_args
step=args.start_from
sl=args.stepslist

steplist=["hifiasm","purge_haplotigs","short_read_filtering","short_read_mapping","pilon","juicer","3DDNA"]
if sl:
    print("hifiasm (already the start of the pipeline so will not have any meaningful effect)\npurge_haplotigs\nshort_read_filtering\nshort_read_mapping\npilon\njuicer\n3DDNA")
    quit()

if step != None and step not in steplist:
    print("ERROR: please input a valid step to start from (to see a list of possible steps use -sl)")
    quit()

if threads==None:
    threads=94

if outfile.endswith("/"):
    outfile=outfile[:-1]

if snail == None:
    snail="lmaximus"

if juicerdir==None:
    juicerdir="/media/data2/juicer-1.6/CPU"

if dnadir==None:
    dnadir="/media/data2/3d-dna-master_altered_fixed"

if pilonjar==None:
    pilonjar="/home/milan/pilon-1.24.jar"

if os.path.isdir(juicerdir) == False:
    print(f"{juicerdir} does not exist, exiting")

if os.path.isdir(dnadir) == False:
    print(f"{dnadir} does not exist, exiting")

if os.path.isfile(pilonjar) == False:
    print(f"{pilonjar} does not exist, exiting")

if hifiasm_args==None:
    hifiasm_args=""
if pilon_args==None:
    pilon_args=""
if juicer_args==None:
    juicer_args=""
if dna3d_args==None:
    dna3d_args=""

if os.path.isdir(outfile) and step==None:
    while True:
        user_input = input(f"directory {outfile} already exist, would you like to overwrite this directory? (yes/no):")
        if user_input.lower() in ["yes", "y"]:
            shutil.rmtree(outfile)
            outpath=outfile
            break
        elif user_input.lower() in ["no", "n"]:
            print("changing output directory")
            val=2
            while True:
                newpath=f"{outfile}_{str(val)}"
                if os.path.isdir(newpath) == False:
                    outpath=newpath
                    break
                else:
                    val+=1
            break
        else:
            print("Invalid input. Please enter yes/no.")
else:
    outpath=outfile

workdir=outpath+"/working_directory"
resultsdir=outpath+"/results"
assemblydir=f"{workdir}/assembly"

if os.path.isdir(outpath) == False:
    os.mkdir(outpath)
    os.mkdir(workdir)
    os.mkdir(resultsdir)
    os.mkdir(assemblydir)



if step==None:
    stepnum=0
else:
    stepnum=steplist.index(step)
deactivate=range(0,stepnum+1)
stepdir={}
for num in range (1,8):
    if num in deactivate:
        stepdir[num]=False
    else:
        stepdir[num]=True

def stepscheck(step):
    if stepdir[step]:
        return True
    else:
        return False






steps=1

if os.path.isdir(f"{assemblydir}/hifiasm")==False:
    os.mkdir(f"{assemblydir}/hifiasm")

if stepscheck(steps):
    if pac:
        #assemble pacbio data
        outpath=f"{assemblydir}/hifiasm/{snail}"
        subprocess.run(f"hifiasm {hifiasm_args} -o {outpath} -t {threads} {longr}", shell=True)
        assemblygraph = f"{outpath}.bp.p_ctg.gfa"
        #convert to fasta
        outpath=f"{assemblydir}/hifiasm/{snail}_assembly.fa"
        command="""awk '/^S/{print ">"$2"\\n"$3}' """+assemblygraph+" | fold > "+outpath
        subprocess.run(command, shell=True)
        assembly=outpath

        #map reads to assembly itself
        outpath=f"{assemblydir}/hifiasm/{snail}_selfalignment.bam"
        subprocess.run(f"pbmm2 align {assembly} {longr} {outpath}",shell=True)
        selfalignment=outpath

        #sort bam output
        outpath=selfalignment.replace(".bam","_sorted.bam")
        subprocess.run(f"samtools sort -@ {threads} {selfalignment} -o {outpath}", shell=True)
        assembly_alignment=outpath

    else:
        #predict genome size using kmerfreq and gce
        #subprocess.run(f"kmerfreq  -k 17 {shortr_cleaned}", shell=True)
        #subprocess.run(f"gce -f {kmerfreqout} -g {kmernum???}", shell=True)
        #size=gceout???

        print("sorry this feature is not supported yet, please use pacbio data")
        quit()

        #filter nanopore data
        outpath= f"{assemblydir}/{snail}_nanopore_filtered"
        #subprocess.run(f"chopper -q 10 -i {longr} > {outpath}", shell=True)
        longr_filtered=outpath

        #assemble nanopore data
        outpath= f"{assemblydir}"
        subprocess.run(f"flye --nano-raw {longr_filtered} --out-dir {outpath}", shell=True)
        assembly= f"{outpath}/"
        #assembly_alignment = ???

    outpath=f"{resultsdir}/{snail}_busco_raw_assembly"
    subprocess.run(f"busco -i {assembly} -f -m genome -c {threads} -l mollusca -f -o {outpath}", shell=True)
else:
    if pac:
        outpath=f"{assemblydir}/hifiasm/{snail}"
        assemblygraph = f"{outpath}.bp.p_ctg.gfa"
        outpath=f"{assemblydir}/hifiasm/{snail}_assembly.fa"
        assembly=outpath
        outpath=f"{assemblydir}/hifiasm/{snail}_selfalignment.bam"
        selfalignment=outpath
        outpath=selfalignment.replace(".bam","_sorted.bam")
        assembly_alignment=outpath
    else:
        print("sorry this feature is not supported yet, please use pacbio data")
        quit()

steps+=1

if os.path.isdir(f"{assemblydir}/purge_haplotigs")==False:
    os.mkdir(f"{assemblydir}/purge_haplotigs")

if stepscheck(steps):
    #remove haplotigs from the assembly
    subprocess.run(f"purge_haplotigs  readhist -b {assembly_alignment} -g {assembly}", shell=True)
    genecov=f"{assembly_alignment.split('/')[-1]}.200.gencov"
    coveragehist=f"{assembly_alignment.split('/')[-1]}.histogram.200.png"
    while True:
        user_input = input(f"please input the low cutoff value from the generated histogram located at {coveragehist}:")
        try:
            int(user_input)
        except:
            print("please input an integer value")
        else:
            low=user_input
            break
    while True:
        user_input = input(f"please input the high cutoff value from the generated histogram located at {coveragehist}:")
        try:
            int(user_input)
        except:
            print("please input an integer value")
        else:
            high=user_input
            break
    while True:
        user_input = input(f"please input the location of the lowest point between the two peaks on the generated histogram located at {coveragehist}:")
        try:
            int(user_input)
        except:
            print("please input an integer value")
        else:
            mid=user_input
            break

    outpath=f"{assemblydir}/purge_haplotigs/{snail}_coverage_stats.csv"
    subprocess.run(f"purge_haplotigs cov -i {genecov} -l {low} -m {mid} -h {high} -o {outpath}", shell=True)
    coverage_stats=outpath

    outpath=f"{assemblydir}/purge_haplotigs/{snail}_purged"
    subprocess.run(f"purge_haplotigs  purge  -g {assembly}  -c {coverage_stats} -o {outpath}", shell=True)
    assembly_purged=f"{outpath}.fasta"
else:
    genecov=f"{assembly_alignment}.200.genecov"
    coveragehist=f"{assembly_alignment}.histogram.200.png"
    outpath=f"{assemblydir}/purge_haplotigs/{snail}_coverage_stats.csv"
    coverage_stats=outpath
    outpath=f"{assemblydir}/purge_haplotigs/{snail}_purged"
    assembly_purged=f"{outpath}.fasta"

if ns==False:
    if os.path.isdir(f"{assemblydir}/clean_illumina")==False:
        os.mkdir(f"{assemblydir}/clean_illumina")

    steps+=1
    if stepscheck(steps):
        #filter short read data (clean_illumina must be manually installed)
        outpath=f"{assemblydir}/clean_illumina/{snail}_short_filtered_1"
        outdata=outpath+"_report.txt"
        outpath=outpath+".fasta"
        subprocess.run(f"clean_lowqual {shortr1} {outpath} {outdata} -t {threads}", shell=True)
        shortr1filt=outpath

        outpath=f"{assemblydir}/clean_illumina/{snail}_short_filtered_2"
        outdata=outpath+"_report.txt"
        outpath=outpath+".fasta"
        subprocess.run(f"clean_lowqual {shortr2} {outpath} {outdata} -t {threads}", shell=True)
        shortr2filt=outpath
    else:
        outpath=f"{assemblydir}/clean_illumina/{snail}_short_filtered_1"
        outdata=outpath+"_report.txt"
        outpath=outpath+".fasta"
        shortr1filt=outpath

        outpath=f"{assemblydir}/clean_illumina/{snail}_short_filtered_2"
        outdata=outpath+"_report.txt"
        outpath=outpath+".fasta"
        shortr2filt=outpath

    if os.path.isdir(f"{assemblydir}/bwa_mem")==False:
        os.mkdir(f"{assemblydir}/bwa_mem")

    steps+=1
    if stepscheck(steps):
        #map short read data assembly
        outpath=f"{assemblydir}/bwa_mem/{snail}_alignment.sam"
        subprocess.run(f"bwa index {assembly_purged}", shell=True)
        subprocess.run(f"bwa mem {assembly_purged} {shortr1filt} {shortr2filt} -t {threads}  > {outpath}" , shell=True)
        alignment=outpath

        #Convert short read data to indexed bam file
        subprocess.run(f"samtools view -@ {threads} -b -S {outpath} > {outpath.replace('.sam','.bam')}", shell=True)
        subprocess.run(f"samtools sort -@ {threads} {outpath.replace('.sam','.bam')} -o {outpath.replace('.sam','_sorted.bam')}", shell=True)
        subprocess.run(f"samtools index -@ {threads} {outpath.replace('.sam','_sorted.bam')}", shell=True)
        alignmentsort=outpath.replace('.sam','_sorted.bam')
        os.remove(alignment)
        os.remove(alignment.replace('.sam','.bam'))
    else:
        outpath=f"{assemblydir}/bwa_mem/{snail}_alignment.sam"
        alignment=outpath
        alignmentsort=outpath.replace('.sam','_sorted.bam')

    if os.path.isdir(f"{assemblydir}/pilon")==False:
        os.mkdir(f"{assemblydir}/pilon")

    #polish assembly (pilon must be manually installed))
    steps+=1
    if stepscheck(steps):
        outpath=f"{assemblydir}/pilon/{snail}_corrected_assembly"
        subprocess.run(f"java -Xmx900G -jar {pilonjar}  --frags {alignmentsort} --genome {assembly_purged} --outdir {assemblydir}/pilon --output {outpath.split('/')[-1]} --fix bases --nonpf --minqual 20 {pilon_args}", shell=True)
        assembly_corrected=outpath+".fasta"
        shutil.copyfile(assembly_corrected, f"{juicerdir}/corrected_assembly.fastq")
        subprocess.run(f"bwa index {juicerdir}/corrected_assembly.fastq", shell=True)
    else:
        outpath=f"{assemblydir}/pilon/{snail}_corrected_assembly"
        assembly_corrected=outpath+".fasta"

else:
    assembly_corrected=assembly_purged

if sj:
    print("-sj flag was triggered, quitting")
    quit()

if os.path.isdir(f"{assemblydir}/juicer")==False:
    os.mkdir(f"{assemblydir}/juicer")

steps+=1
if stepscheck(steps):
#map hi-c data to assembly (Juicer must be manually installed)
    if os.path.isdir(f"{juicerdir}/aligned"):
        shutil.rmtree(f"{juicerdir}/aligned")
    if os.path.isdir(f"{juicerdir}/splits"):
        shutil.rmtree(f"{juicerdir}/splits")
    os.mkdir(f"{juicerdir}/splits")
    shutil.copyfile(hicr1, f"{juicerdir}/splits/{snail}_HIC_R1.fastq.gz")
    shutil.copyfile(hicr2, f"{juicerdir}/splits/{snail}_HIC_R2.fastq.gz")
    subprocess.run(f"bash juicer.sh -z corrected_assembly.fastq {juicer_args} -p chrom.sizes -S early", cwd=juicerdir, shell=True)
    shutil.rmtree(f"{juicerdir}/splits")
    shutil.copyfile(f"{juicerdir}/aligned/merged_nodups.txt",f"{dnadir}/merged_nodups.txt")
    shutil.copyfile(f"{juicerdir}/aligned/merged_nodups.txt",f"{assemblydir}/{juicer}/{snail}_merged_nodups.txt")

steps+=1
if stepscheck(steps):
    #create hi-c map and final assembly (3DDNA must be manually installed)
    shutil.copyfile(assembly_corrected, f"{dnadir}/{assembly_corrected.split('/')[-1]}")
    subprocess.run(f"bash run-asm-pipeline.sh -r 0 {dna3d_args} {assembly_corrected.split('/')[-1]} merged_nodups.txt", cwd=dnadir, shell=True)
    finalassembly=f"""{dnadir}/{assembly_corrected.split("/")[-1].replace(".fasta",".FINAL.fasta")}"""
    finalassemblyasm=f"""{dnadir}/{assembly_corrected.split("/")[-1].replace(".fasta",".FINAL.assembly")}"""
    finalassemblyhic=f"""{dnadir}/{assembly_corrected.split("/")[-1].replace(".fasta",".final.hic")}"""
    shutil.copyfile(finalassembly, f"{resultsdir}/{snail}_final_assembly_uncorrected.fasta")
    shutil.copyfile(finalassemblyasm, f"{resultsdir}/{snail}_final_assembly_uncorrected.assembly")
    shutil.copyfile(finalassemblyhic, f"{resultsdir}/{snail}_final_assembly_uncorrected.hic")
    outpath=f"{resultsdir}/{snail}_busco_final_assembly"
    subprocess.run(f"busco -i {finalassembly}  -m genome -f -c {threads} -l mollusca -o {outpath}", shell=True)

print("assembly complete, now check your hi-c heatmap and rerun 3DDNA with different parameters if required, after achieving a satisfactory heatmap manually fix any remainging errors using juicebox and place this corrected assembly in the results folder with the name (snail)_final_assembly_corrected.fasta")






