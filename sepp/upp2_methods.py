import os, glob
import subprocess

from sepp.filemgr import get_root_temp_dir
from sepp.hmm_concurrent import *
from sepp.scheduler import JobPool

def create_dirs(path_name):
    if not os.path.exists(path_name): 
        subprocess.call(['mkdir', path_name])

def parse_strats(stratsfile): 
    stratlist = []
    with open(stratsfile, "r") as reader:
        lines = reader.readlines()
        for line in lines:
            stratlist.append(line.strip())
    return stratlist


def makedirstruct(dirpath): 
    print("[makedirstruct]")
    # clear the folders
    dirpath = os.path.join(dirpath, "ensembleData")
    create_dirs(dirpath)
    for firstlvl in ["trueAlignment","temporaryFileSave","subsetTrueAln", "Searcher","queryToHmm","newHMM","ML","initialHMM","hmmSeqAlign","hmmScores","hmmQueryList","fullPrediction", "seqFileNames"]:
        
        firstdir = '%s/%s' % (dirpath, firstlvl)
        if os.path.isdir(firstdir):
            subprocess.call(['rm', '-rf', firstdir])

        subprocess.call(['mkdir', firstdir])
        if firstlvl == 'fullPrediction':
            subprocess.call(['mkdir', firstdir + '/sequences/'])
        elif firstlvl == 'hmmQueryList':
            for z in ['inputQuery', 'merged', 'predictedQuery']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'hmmScores':
            for z in ['newRawOld', 'processedOld', 'rawOld', 'scoresFull', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'ML':
            for z in ['hmm', 'models', 'queryHMM', 'scores', 'sequences', 'temporaryStorage']:
                subprocess.call(['mkdir', firstdir + '/' + z])
            subprocess.call(['mkdir', firstdir + '/sequences/np'])
        elif firstlvl == 'newHMM':
            for z in ['columnSets', 'hmm', 'newHMMseq']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'queryToHmm':
            for z in ['original', 'withUPP']:
                subprocess.call(['mkdir', firstdir + '/' + z])
        elif firstlvl == 'Searcher':
            subprocess.call(['mkdir', '%s/scoreFiles' % firstdir])
        elif firstlvl == 'trueAlignment':
            for z in ['original', 'subset']:
                subprocess.call(['mkdir', firstdir + '/' + z])

def run_upp_strats(abstract_algorithm, dirname, hier_upp, adjusted_bitscore, doResort=False):
    ''' Note: abstract_algorithm can be used to do commands like below.
        The moment you call addHMMBuildJob, jobs will get enqueued and run eventually
    addHMMBuildJob(abstract_algorithm, <hmmbuild profile output file path>, <fasta file path>)
    JobPool().wait_for_all_jobs()
    '''
    print("[run_upp_strats]")

    # globdir = glob.glob(os.path.join(dirname, "output*"))
    # assert len(globdir) == 1
    # outputdirname = globdir[0]
    root_temp_dir = get_root_temp_dir()
    hmm_folder_prefix = "%s/root/P_0/" % root_temp_dir
    # Nit: We should be able to retrieve the same information from options().fragment_file
    fragment_file_suffix = os.listdir(root_temp_dir + "/fragment_chunks/")[0]
    fragment_file = "%s/fragment_chunks/%s" % (root_temp_dir, fragment_file_suffix)
    # predictionName = './%s/UPPoutput/%s_output_alignment.fasta' % dirpath

    # TODO: Prediction name may get passed in
    prediction_name = ""
    # TODO: this is something we need to remove
    true_alignment = ""
    dataset_name = "default-value-not-empty"

    set_all_file_names(hmm_folder_prefix, fragment_file, true_alignment, prediction_name, root_temp_dir, dataset_name)
    save_initial_steps(abstract_algorithm)
    if hier_upp:
        strat = 'stefan_fastUPP'
        hierchySearch(abstract_algorithm)
    elif adjusted_bitscore:
        strat = 'stefan_UPPadjusted'

    print("[processing %s]" % strat)
    print("[running scoresToHMMSeq]")
    scoresToHMMSeq(strat)
    print("[running buildAlignMerge, doResort is %s]" % doResort)
    buildAlignMerge(abstract_algorithm, strat, doResort=doResort)

