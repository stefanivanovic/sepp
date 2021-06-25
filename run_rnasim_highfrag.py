import os, glob
import subprocess

indata = "/home/gchu4/scratch/warnow/ensembles/dataset/UnalignFragTree/high_frag/RNASim/" 
root_tmpdir = "/home/gchu4/scratch/warnow/ensembles/unweighted_output/RNASimHF"
packagedir = "/home/gchu4/scratch/warnow/ensembles/sepp/tools/bundled/Linux/"

for r in range(12, 20):    
    rep = "R"+str(r)
    unaligned_frag_dir = indata + rep + "/unaligned_frag.txt"
    unaligned_frag = glob.glob(unaligned_frag_dir)[0] 
    assert len(glob.glob(unaligned_frag_dir)) == 1
    print("Running...", unaligned_frag)

    for decomp_size in [10, 2]: 
        print("Processing", rep + "_decomp" + str(decomp_size))
        currdir = indata + rep + "/"
        if not os.path.exists(root_tmpdir + "/" + rep):
            subprocess.call(['mkdir', root_tmpdir + "/" + rep])

        tmpdir = root_tmpdir + "/" + rep + "/" + str(decomp_size)

        alignfile = currdir + "pasta_backbone_align.txt"
        treefile = currdir + "pasta_backbone.tre"

        if not os.path.exists(tmpdir):
            subprocess.call(['mkdir', tmpdir])

        configfile = root_tmpdir + "/" + rep + "/" + "decomp" + str(decomp_size) + ".config"
        with open(configfile, "w+") as f:
            f.write('[commandline]\n')
            f.write("alignment=%s\n" % alignfile)
            f.write("tree=%s\n" % treefile)
            f.write("sequence_file=%s\n" % unaligned_frag)
            f.write("alignmentSize=%s\n" % str(decomp_size))
            f.write("molecule=dna\n")
            f.write("tempdir=%s\n" % tmpdir)
            f.write("cpu=1\n")
            f.write('[upp2]\n')
            f.write('decomp_only=False\n')
            f.write('bitscore_adjust=False\n')
            f.write('hier_upp=False\n')
            f.write('early_stop=False\n')
            f.write('[hmmsearch]\n')
            f.write('user_options= --max --tblout \n')
            #f.write("outdir=output\n")

        #subprocess.call(['python', 'run_upp.py', '-s', unaligned_frag, '-x', '1', '-a', alignfile, '-t', treefile, '-A', str(decomp_size), '-c', configfile, '--tempdir', tmpdir]) # -B 200

        subprocess.call(['python', 'run_upp.py', "-c", configfile])

        outputfile = root_tmpdir + "/" + rep + "decomp" + str(decomp_size)
        subprocess.call(['mv', 'output_alignment.fasta', outputfile + "_output_alignment" + ".fasta"])
        subprocess.call(['mv', 'output_alignment_masked.fasta', outputfile + "_output_alignment_masked.fasta"])
        subprocess.call(['mv', 'output_insertion_columns.txt', outputfile + "_output_insertion_columns.txt"])

