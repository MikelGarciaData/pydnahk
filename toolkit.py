import itertools
import numpy as np
from pydna.genbank import Genbank
from Bio import SeqIO
from Bio.Restriction import AllEnzymes
from Bio.Restriction import AllEnzymes

def feature_table(sequence):
    sorted_feat = sorted(sequence.features, key=lambda x: int(x.location.start))
    
    # misc_feature is MCS?
    print(f"{'feature':15} {'start':6}{'end':6} {'notes'}")
    for feature in sorted_feat:
        start = feature.location.start
        end = feature.location.end
        print(f"{feature.type:15} {start:6d}{end:6d} {feature.qualifiers.get("note")}")
        #print(feature.type, feature.location, feature.qualifiers.get("note"))


def load_from_genebank(email, sequence):
    gb = Genbank(users_email=email)  
    seq = gb.nucleotide(sequence)   
    return seq


def cut_enzyme_info(seq, enzyme_list=AllEnzymes):
    cut_sites = {}
    for enzyme in enzyme_list:
        cut = enzyme.search(seq)
        if cut:
            cut_sites[enzyme] = cut
    return cut_sites


def _get_mcs_cuts(sequence, mcs_feat):
    # Get location of MCS
    locs = mcs_feat.location
    start, end = locs.start, locs.end 
    mcs = sequence.seq[start: end]
    
    # Find enzyme cuts
    enz_cuts = cut_enzyme_info(mcs)

    # Get enzymes in the MCS
    mcs_enz = list(enz_cuts.keys())
    
    # Get cuts and add offset
    mcs_cuts = np.array(list(set(itertools.chain.from_iterable(list(enz_cuts.values()))))) + start
    return {"enzymes": mcs_enz, "cuts": mcs_cuts}


def get_mcs_cuts(sequence):
    mcs_list = list(filter(lambda feat: feat.type=='misc_feature', sequence.features))
    mcs_enz_cut = []
    for mcs in mcs_list:
        mcs_enz_cut.append(_get_mcs_cuts(sequence, mcs))
    return mcs_enz_cut



def get_non_mcs_regions(sequence):
    # get non-mcs regions
    mcs_list = list(filter(lambda feat: feat.type=='misc_feature', sequence.features))
    mcs_list
    sequence = sequence
    landmarks = []
    for mcs in mcs_list: 
        locs = mcs.location
        start, end = locs.start, locs.end
        landmarks.append(start)
        landmarks.append(end) 
    
    landmarks = np.array(sorted(landmarks))
    landmarks = np.roll(landmarks, shift=1)
    landmarks
    
    non_mcs_regions = []
    for i in range(0, len(landmarks), 2):
        start, end = landmarks[i:i+2]
        non_mcs = sequence.seq[start: end]
        non_mcs_regions.append(non_mcs)

    return non_mcs_regions

