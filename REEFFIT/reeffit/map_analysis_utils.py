from collections import defaultdict
# import pdb

from matplotlib.pylab import *
from scipy.stats import stats
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as hcluster

import rdatkit.secstr as ss

from reactivity_distributions import *

import joblib

k = 0.0019872041
T = 310.15 # temperature


def _AICc(lhood, nparams, ndata):
    """Akaike information criterion (corrected form) for fit assessment and (in some cases) model selection."""
    return 2 * nparams - 2 * lhood + (2 * nparams * (nparams + 1)) / (ndata - nparams -1)


def AICc(lhood, data, W, E_d, E_c, Psi):
    """
    Akaike information criterion for the REEFFIT factor mode.

    Args:
        lhood (float): The likelihood value from the model
        data (numpy array): The chemical mapping data
        W (numpy array): The structure weights
        E_d (numpy array): The structure latent reactivities
        E_c (numpy array): The structure latent contact matrices
        Psi (numpy array): The sequence-dependent covariance matrices

    Returns
    """
    return _AICc(lhood, W.size + E_d.size + Psi.shape[0] + E_c[logical_not(isnan(E_c))].size, data.size)

def _acc(typ1, typ2, cstart=0, cend=-1):
    typ1a=typ1.replace("(","+")
    typ1a=typ1a.replace(")","+")
    typ2a=typ2.replace("(","+")
    typ2a=typ2a.replace(")","+")
    res = 0
    if(cend==-1):
        ce=len(typ1)
    else:
        ce=cend
    for i in xrange(cstart,ce):
        if(typ1a[i]!=typ2a[i]):
            res+=1
    return(res)


# Calculate mutual information between types
def _mutinf(typ1, typ2, cstart=0, cend=-1):
    #TODO: Right now (just using paired and unpaired distributions)
    # we just need to check if the two types are the same or not
    # in the future, when we have confident distributions between all structure subtypes
    # (bulges, internal pairs, etc.) we will replace this with the true mutual information
    # between types -- given that these structure subtypes are fixed, we could even cache them!
    res = 0
    if(cend==-1):
        ce=len(typ1)
    else:
        ce=cend
#    for i in xrange(len(typ1)):
    for i in xrange(cstart, ce):
        if typ1[i] != typ2[i]:
            res += 1
    return res


def bpdist(sbp1, sbp2):
    res = 0
    for bp1_1, bp1_2 in sbp1:
        for bp2_1, bp2_2 in sbp2:
            if not ((bp1_1 == bp2_1 and bp1_2 == bp2_2) or (bp1_2 == bp2_1 and bp1_1 == bp2_2)):
                res += 1
    return res

def get_struct_types(structures, cutoff=-1):
    struct_types = []
    for i in xrange(len(structures[0])):
        # Just mark paried as p and unpaired as u, later we may refine the model
        struct_types.append(['u' if s[i] == '.' else 'p' for s in structures])
    return struct_types


def get_structure_distance_matrix(structures, struct_types, distance='mutinf', cstart=0, cend=-1):
    nstructs = len(struct_types[0])
    if distance == 'basepair':
        if(cend==-1):
            ce = len(structures[0])
        else:
            ce = cend
        structures_bp = [ss.SecondaryStructure(dbn=struct[cstart:cend]).base_pairs() for struct in structures]
    D = zeros([nstructs, nstructs])
    for i in xrange(nstructs):
        st1 = [s[i] for s in struct_types]
        D[i, i] = 0
        for j in xrange(i+1, nstructs):
            if distance == 'mutinf':
                D[i, j] = _mutinf(st1, [s[j] for s in struct_types], cstart, cend)
            elif distance == 'basepair':
                D[i, j] = bpdist(structures_bp[i], structures_bp[j])
            elif distance == 'acc':
                D[i, j] = _acc(structures[i], structures[j], cstart, cend)
            else:
                raise ValueError('Distance %s not recognized: options are "mutinf" and "basepair"' % distance)
            D[j, i] = D[i, j]
    return D

def get_structure_distance_matrix_2(structures, struct_types, njobs):
    nstructs = len(struct_types[0])
    structures_bp = [ss.SecondaryStructure(dbn=struct).base_pairs() for struct in structures]
    D = zeros([nstructs, nstructs])

    def dist_fun(i):
#        print 'Calculating dist for row %s' % i
        dd = zeros([1, nstructs])
        for j in xrange(i+1, nstructs):
            dd[0, j] = bpdist(structures_bp[i], structures_bp[j])
        return dd
    sys.modules[__name__].dist_fun = dist_fun
    res = joblib.Parallel(n_jobs=njobs)(joblib.delayed(dist_fun)(i) for i in xrange(nstructs))
    for k in xrange(nstructs):
        D[k,:] = res[k]
    D = D + D.T - diag(diag(D))
    print 'Dist calculation done'
    return D


#def get_structure_distance_matrix_2(structures, struct_types):
#    cdef struct pair_t:
#        int pos_x
#        int pos_y
#    cdef void make_pairs(list pairs, pair_t *cpairs, int num_pairs)

#    nstructs = len(struct_types[0])
#    structures_bp_mat = zeros([nstructs, len(structures[0]), len(structures[0])])
#    for k in xrange(nstructs):
#        bps = ss.SecondaryStructure(dbn=structures[k]).base_pairs()
#        for bp in bps:
#            structures_bp_mat[k, bp[0], bp[1]] = 1
#            structures_bp_mat[k, bp[1], bp[0]] = 1
#    
#    D = zeros([nstructs, nstructs])
#
#    for i in xrange(nstructs):
#        D[i, i] = 0
#        for j in xrange(i+1, nstructs):
#            D[i, j] = sum(absolute(structures_bp_mat[i,:,:] - structures_bp_mat[j,:,:]))
#            D[j, i] = D[i, j]
#        print i, j
 #   return D

def cluster_structures(struct_types, structures=[], distance='mutinf', expected_medoids=-1, cstart=0, cend=-1):
    nstructs = len(struct_types[0])
    print 'Clustering structures'
    print expected_medoids
    # Ok, find a subset of structures that are mostly "orthogonal", clustering
    # with mutual information (mutinf option) or base pair distance (basepair option)
    # Looks ugly, but seems that numpy's symmetric matrix support is limited (i.e. none)
    print (cstart, cend)
    D = get_structure_distance_matrix(structures, struct_types, distance='acc', cstart=cstart, cend=cend)
    D3 = get_structure_distance_matrix(structures, struct_types, distance='mutinf', cstart=cstart, cend=cend)
    D2 = get_structure_distance_matrix(structures, struct_types, distance='basepair', cstart=cstart, cend=cend)
    # Proxy for ensemble centroid (EC)
    minsumdist = inf
    for i in xrange(nstructs):
        sumdist = D[i, :].sum()
        if sumdist < minsumdist:
            minsumdist = sumdist
            EC = i
    Z = hcluster.linkage(squareform(D), method='complete')

    prevnumclusts = nstructs
    maxCH = -inf
    assignments = dict([(ns, [ns]) for ns in xrange(nstructs)])
    maxmedoids = dict([(ns, ns) for ns in xrange(nstructs)])
    for t in arange(0, D.max(), D.max()/10.):
        flatclusters = hcluster.fcluster(Z, t, criterion='distance')
        clusters = defaultdict(list)
        for elem, cidx in enumerate(flatclusters):
            clusters[cidx].append(elem)
        numclusts = len(clusters)
        if numclusts == prevnumclusts:
            continue
        else:
            prevnumclusts = numclusts
#        if numclusts == 1:
        if numclusts <= expected_medoids:
            break
        # Choose cluster medoids
        medoids = {}
        mindist = 0
        for cidx, elems in clusters.iteritems():
            for j, elem in enumerate(elems):
                sumdist = D[elem, elems].sum() / 2.
                if j == 0 or sumdist < mindist:
                    mindist = sumdist
                    medoids[cidx] = elem

        #if expected_medoids > 0:
            # If we already estimated the expected structures,
            # then just retrieve the cluster assignments and best medoids
#            maxmedoids = medoids
#            assignments = clusters
#            break

        # Calculate Calinki-Harabasz index to choose clusters
        # Using between-cluster distances (BC) and within-cluster (WC) distances
        BC = sum([len(elems) * D[medoids[cidx], EC] for cidx, elems in clusters.iteritems() if len(elems) > 1])
        WC = sum([D[medoids[cidx], elems].sum() / 2. for cidx, elems in clusters.iteritems() if len(elems) > 1])
        if BC == 0 or WC == 0:
            CH = -inf
        else:
            CH = (BC / (numclusts - 1)) / (WC / (D.shape[0] - numclusts))
        if CH > maxCH:
            maxCH = CH
            maxmedoids = medoids
            assignments = clusters
        print 'For %s threshold we have %s clusters with %s CH' % (t, numclusts, CH)
    print 'Done clustering structures, with %s clusters' % numclusts

#    print assignments
    return maxmedoids, assignments, D, D2


def normalize(data):
    wtdata = array(data)
    wtdata[wtdata < 0] = 0
    q1 = stats.scoreatpercentile(wtdata, 25)
    q3 = stats.scoreatpercentile(wtdata, 75)
    interquart = q3 - q1
    tenperc = stats.scoreatpercentile(wtdata[wtdata <= q1 + 1.5 * interquart], 90)
    maxav = wtdata[wtdata >= tenperc].mean()
    wtdata = wtdata / maxav
    return wtdata


# Energies is the energies used to calculate the weights,
# clusters is a dictionary with lists of structure indices
def calculate_weights(energies, clusters=None):
    if clusters is not None:
        W = zeros([energies.shape[0], len(clusters)])
        struct_weights_by_clust = {}
        part_func = zeros([energies.shape[0]])
        for c, struct_indices in clusters.iteritems():
            struct_weights_by_clust[c] = zeros([energies.shape[0], len(struct_indices)])
            for i, j in enumerate(struct_indices):
                struct_weights_by_clust[c][:, i] = exp(-energies[:, j] / (k * T))
                part_func[:] += struct_weights_by_clust[c][:, i]
        for i, c in enumerate(struct_weights_by_clust):
            for j in xrange(energies.shape[0]):
                struct_weights_by_clust[c][j, :] /= part_func[j]
            W[:, i] = struct_weights_by_clust[c].sum(axis=1)
        return W, struct_weights_by_clust
    else:
        W = zeros(energies.shape)
        for j in xrange(energies.shape[0]):
            W[j, :] = exp(-energies[j, :] / (k * T)) / exp(-energies[j, :] / (k * T)).sum()
        return W


def mock_data(sequences, structures=None, energy_mu=0.5, energy_sigma=0.5, obs_sigma=0.01, paired_sampling=lambda : SHAPE_paired_sample(), unpaired_sampling= lambda : SHAPE_unpaired_sample(), contact_sampling=lambda : SHAPE_contacts_diff_sample(), mutpos=None, c_size=3, return_steps=False, correlate_regions=False):
    if structures is not None:
        print 'Generating mock data'
        print 'Getting "true" free energies (from RNAstructure)'
        true_energies = get_free_energy_matrix(structures, sequences)

        print 'Energies'
        print true_energies

        # Add some noise to the energies
        noise = energy_sigma*randn(true_energies.shape[0], true_energies.shape[1]) + energy_mu
        #noise = rand(true_energies.shape[0], true_energies.shape[1])
        noised_energies = true_energies + noise
        weights_noised = calculate_weights(noised_energies)
        weights = calculate_weights(true_energies)
        print 'Weights'
        print weights
        print 'Noised Weights'
        print weights_noised

        print 'Generating mock reactivities'
        # Mock reactivities
        reacts = zeros([len(structures), len(sequences[0])])
        prev_s = ''
        MIN_REACT = 1e-5
        P_CONTACT = 1
        if correlate_regions:
            for j, s in enumerate(structures):
                for i in xrange(len(sequences[0])):
                    if s[i] != prev_s:
                        if s[i] == '.':
                            curr_val = unpaired_sampling()
                        else:
                            curr_val = max(MIN_REACT, paired_sampling())
                            curr_val = paired_sampling()
                    reacts[j, i] = max(MIN_REACT, curr_val + 0.01 * SHAPE_contacts_diff_sample())
                    #reacts[j, i] = curr_val
                    prev_s = s[i]
        else:
            prevstate = '.'
            for j, s in enumerate(structures):
                for i in xrange(len(sequences[0])):
                    if s[i] == '.':
                        reacts[j, i] = unpaired_sampling()
                        prevstate = '.'
                    else:
                        if prevstate == '.' or (i < len(sequences[0]) - 1 and s[i + 1] == '.'):
                            reacts[j, i] = paired_sampling() * 1.5
                        else:
                            reacts[j, i] = paired_sampling() * 0.3
                        prevstate = '('

        data = dot(weights_noised, reacts)
        data_orig = dot(weights_noised, reacts)
        if mutpos:
            """
            for i, pos in enumerate(mutpos):
                if i >=0:
                    for k in xrange(-c_size/2, c_size/2+1):
                        if pos+ k < data.shape[1] and rand() > P_CONTACT:
                            data[i,pos+ k] = contact_sampling()
            """
            print 'Simulate diagonal and off-diagonal contact sites'
            max_tries = 1000
            def add_local_perturb(reactivity, weight):
                dd = weight*contact_sampling()
                tries = 0
                while reactivity + dd < MIN_REACT or reactivity + dd > 4.5:
                    dd = weight*contact_sampling() * 0.1
                    tries += 1
                    if tries > max_tries:
                        if reactivity + dd > 4.5:
                            return 4.5
                        if reactivity + dd < MIN_REACT:
                            return MIN_REACT - reactivity
                return dd

            bp_dicts = []
            for s, struct in enumerate(structures):
                bp_dicts.append(ss.SecondaryStructure(dbn=struct).base_pair_dict())
            for j in xrange(data.shape[0]):
                for k in xrange(-(c_size - 1) / 2, (c_size - 1) / 2 + 1):
                    for s in xrange(len(structures)):
                        #if weights_noised[j, s] < 0.1:
                        #    continue
                        if type(mutpos[j]) is list:
                            for m in mutpos[j]:
                                if m + k < data.shape[1] and rand() < P_CONTACT:
                                    if structures[s][m + k] == '.':
                                        dd = 0.1*add_local_perturb(data_orig[j, m + k], weights_noised[j, s])
                                    else:
                                        dd = 0.3*add_local_perturb(data_orig[j, m + k], weights_noised[j, s])
                                    if k != 0:
                                        dd *= 0.2
                                    data[j, m + k] += dd
                                if m in bp_dicts[s] and bp_dicts[s][m] + k < data.shape[1] and rand() < P_CONTACT:
                                    data[j, bp_dicts[s][m] + k] += dd

                        else:
                            if mutpos[j] + k < data.shape[1] and rand() < P_CONTACT:
                                if structures[s][mutpos[j]+ k] == '.':
                                    dd = 0.1 * add_local_perturb(data_orig[j, mutpos[j] + k], weights_noised[j, s])
                                else:
                                    dd = 0.3 * add_local_perturb(data_orig[j, mutpos[j] + k], weights_noised[j, s])
                                if k != 0:
                                    dd *= 0.2
                                data[j, mutpos[j] + k] += dd
                            if mutpos[j] in bp_dicts[s] and bp_dicts[s][mutpos[j]] + k < data.shape[1] and rand() < P_CONTACT:
                                data[j, bp_dicts[s][mutpos[j]] + k] += dd


        print 'Adding observational noise'
        data_noised = zeros(data.shape)
        obs_noise_sigmas = []
        params = (0.10524313598815455, 0.034741986764665007)
        for i in xrange(data.shape[1]):
            sigma = rand() * 0.2 + obs_sigma
            sigma = max(stats.distributions.cauchy(loc=params[0], scale=params[1]).rvs(), 0.001) * 0.2
            data_noised[:, i] = data[:, i] + randn(data.shape[0]) * sigma
            obs_noise_sigmas.append(sigma)
        data_noised = data + randn(data.shape[0], data.shape[1]) * obs_sigma
        if return_steps:
            return dot(weights, reacts), data, data_noised, true_energies, weights_noised, reacts, obs_noise_sigmas
        else:
            return data_noised
    else:
        data = zeros(len(sequences), len(sequences[0]))
        for j, seq in enumerate(sequences):
            bppm = ss.partition(sequence)
            unpaired_probs = 1 - bppm.sum(axis=0)
            for i, up in enumerate(unpaired_probs):
                data[j, i] = up*unpaired_sampling() + (1 - up) * paired_sampling()
        data_noised = data + obs_sigma * randn(data.shape[0], data.shape[1])
        return data_noised


def remove_non_cannonical(structure, sequence):
    cannonical_bp = [('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G'), ('A', 'U'), ('U', 'A')]
    bp_dict = ss.SecondaryStructure(dbn=structure).base_pair_dict()
    res_struct = ['.']*len(sequence)
    for n1, n2 in bp_dict.iteritems():
        if (sequence[n1], sequence[n2]) in cannonical_bp:
            if n1 < n2:
                res_struct[n1] = '('
                res_struct[n2] = ')'
            else:
                res_struct[n1] = ')'
                res_struct[n2] = '('
    return ''.join(res_struct)


def get_free_energy_matrix(structures, sequences, algorithm='rnastructure', njobs=1):
    energies = zeros([len(sequences), len(structures)])
    if algorithm == 'uniform':
        energies += 1
    else:

        def efn_fun(seq):
            print 'Calculating structure energies for sequence %s' % seq
            energy = zeros([1, len(structures)])
            energy[0, :] = array(ss.get_structure_energies(seq, [ss.SecondaryStructure(dbn=remove_non_cannonical(s, seq)) for s in structures], algorithm=algorithm))
#            minenergy = energy[1, :].min()
            for i in xrange(len(structures)):
                energy[0, i] = min(energy[0, i], 200)
            return energy
            
        sys.modules[__name__].efn_fun = efn_fun
        print 'Calculating structure energies'
        res = joblib.Parallel(n_jobs=njobs)(joblib.delayed(efn_fun)(seq) for j, seq in enumerate(sequences))
        for j in xrange(len(sequences)):
            energies[j,:] = res[j]

#        for j, seq in enumerate(sequences):
#            print 'Calculating structure energies for sequence %s: %s' % (j, seq)
#            energies[j, :] = array(ss.get_structure_energies(seq, [ss.SecondaryStructure(dbn=remove_non_cannonical(s, seq)) for s in structures], algorithm=algorithm))
#            minenergy = energies[j, :].min()
#            #energies[j, :] -= minenergy
#            for i in xrange(len(structures)):
##                energies[j, i] = max(energies[j, i], 200)
#                energies[j, i] = min(energies[j, i], 200)

##            print energies[j,:]
    return energies


def get_contact_sites(structures, mutpos, nmeas, npos, c_size, restrict_range=None):
    bp_dicts = []
    nstructs = len(structures)
    if restrict_range:
        mutpos_cutoff = [[m + restrict_range[0] if m > 0 else m for m in pos] for pos in mutpos]
    else:
        mutpos_cutoff = mutpos
    for s, struct in enumerate(structures):
        bp_dicts.append(ss.SecondaryStructure(dbn=struct).base_pair_dict())
    contact_sites = {}
    for s in xrange(nstructs):
        contact_sites[s] = zeros([nmeas, npos])
    nstructs = len(structures)
    for j in xrange(nmeas):
        if len(mutpos_cutoff[j]) > 0:
            for m in mutpos_cutoff[j]:
                for s in xrange(nstructs):
                    for k in xrange(-(c_size-1)/2, (c_size-1)/2+1):
                        if m + k >= 0 and m + k < npos:
                            contact_sites[s][j, m + k] = 1
                        if m in bp_dicts[s] and bp_dicts[s][m] + k < npos and bp_dicts[s][m] + k >= 0:
                            contact_sites[s][j, bp_dicts[s][m] + k] = 1
    if restrict_range is not None:
        for s in xrange(nstructs):
            contact_sites[s] = contact_sites[s][:, restrict_range[0]:restrict_range[1]]
    return contact_sites


def get_minimal_overlapping_motif_decomposition(structures, bytype=False, offset=0):
    if type(structures[0]) is str:
        struct_objs = [ss.SecondaryStructure(dbn=s) for s in structures]
    else:
        struct_objs = structures

    def get_motif_id(k, ntlist, pos):
        if bytype:
            return '%s_%s' % (k, pos)
        return '%s_%s' % (k, ';'.join([str(x) for x in ntlist]))

    def get_type_and_ntlist(id):
        typ, ntliststr = id.split('_')
        return typ, [int(x) for x in ntliststr.split(';')]

    pos_motif_map = {}

    elems = [s.explode() for s in struct_objs]
    cover_matrix = ones([len(struct_objs), len(struct_objs[0])])
    motif_ids = []
    for i, s1 in enumerate(struct_objs):
        cover_vec = [False] * len(s1)
        for k, v in elems[i].iteritems():
            for ntlist in v:
                for pos in ntlist:
                    cover_vec[pos] = True
        # Single stranded regions that were not covered by a motif
        # are collapsed into a "motif" we call sstrand
        ssprev = -1
        currssmotif = []
        elems[i]['sstrand'] = []
        foundssmotif = False
        for j in xrange(len(s1)):
            if not cover_vec[j]:
                if ssprev == j - 1:
                    currssmotif.append(j)
                else:
                    if ssprev >= 0:
                        elems[i]['sstrand'].append(currssmotif)
                    currssmotif = [j]
                    foundssmotif = True
                ssprev = j
        if foundssmotif:
            elems[i]['sstrand'].append(currssmotif)

        for k, v in elems[i].iteritems():
            for ntlist in v:
                for pos in ntlist:
                    try:
                        m_idx = motif_ids.index(get_motif_id(k, ntlist, pos+offset))
                        if (pos+offset, m_idx) not in pos_motif_map:
                            pos_motif_map[(pos+offset, m_idx)] = [i]
                        else:
                            if i not in pos_motif_map[(pos + offset, m_idx)]:
                                pos_motif_map[(pos + offset, m_idx)].append(i)
                    except ValueError:
                        motif_ids.append(get_motif_id(k, ntlist, pos+offset))
                        pos_motif_map[(pos + offset, len(motif_ids) - 1)] = [i]

    motif_dist = zeros([len(motif_ids), len(motif_ids)])
    if bytype:
        MAX_DIST = 1.
        MIN_DIST = 1e-5
        for i, mid1 in enumerate(motif_ids):
            for j, mid2 in enumerate(motif_ids):
                typ1, ntlist1 = get_type_and_ntlist(mid1)
                typ2, ntlist2 = get_type_and_ntlist(mid2)
                if typ1 == typ2:
                    nposoverlap = 0.
                    for nt1 in ntlist1:
                        if nt1 in ntlist2:
                            nposoverlap += 1.
                    motif_dist[i, j] = max(MIN_DIST, nposoverlap/max(len(ntlist1), len(ntlist2)))
                else:
                    motif_dist[i, j] = MAX_DIST
                motif_dist[j, i] = motif_dist[i, j]


    """
    for k, d in motif_pos_map.iteritems():
        for i, s1 in enumerate(struct_objs):
            for pos in xrange(len(s1)):
                if pos in d and i in d[pos]:
                    pos_motif_map[pos][k].append(i)

    return motif_pos_map, pos_motif_map, motif_ids
    """
    for i in xrange(cover_matrix.shape[0]):
        for j in xrange(cover_matrix.shape[1]):
            for pos, mid in pos_motif_map.keys():
                if pos == j and i in pos_motif_map[(pos, mid)]:
                    cover_matrix[i, j] = 0
                    break

    return pos_motif_map, motif_ids, motif_dist

    """
    # This is an unfinished implementation of the greedy algorithm for
    # "set covering" applied to minimal motif decomposition
    base_types = []
    bp_dicts = [s.base_pair_dict() for s in struct_objs]
    for i, s in enumerate(struct_objs):
        type_dict = s.explode()
        base_types.append(['']*len(s.dbn))
        for typ, elems in type_dict.iteritems():
            for ntlist in elems:
                for j in ntlist:
                    # Do not discriminate between n-way junctions
                    if 'junction' in typ:
                        base_types[i][j] = 'junction'
                    else:
                        base_types[i][j] = typ

    for i, s in enumerate(struct_objs):
        motif = []
        share_seqs = []
        idx = 0
        while idx < len(s)
            for j in xrange(idx, len(s.dbn)):
                if base_types[i][j] == 'helix':
                    for k, s2 in enumerate(struct_objs):
                        if base_dicts[i][j] == base_dicts[k][j] and k not in share_seqs:
                            if k, not in already_added:
                                share_seqs.append(k)
                                already_added.append(k)
                        if base_dicts[i][j] != base_dicts[k][j] and k in share_seqs:
                            stop_flag = True
                else:
                    for k, s2 in enumerate(struct_objs):
                        if base_types[i][j] == base_types[k][j] and k not in share_seqs:
                            if k not in already_added:
                                share_seqs.append(k)
                        if base_types[i][j] != base_types[k][j] and k in share_seqs:
                            stop_flag = True
                            break
                if stop_flag:
                    break
                else:
                    motif.append(j)
            idx = j
    """


def bpp_matrix_from_structures(structures, weights, weight_err=None, signal_to_noise_cutoff=0, flip=False, symmetric=True):
    npos = len(structures[0])
    bppm = zeros([npos, npos])
    if weight_err is not None:
        bppm_err = zeros([npos, npos])

    for i, s in enumerate(structures):
        for n1, n2 in ss.SecondaryStructure(dbn=s).base_pairs():
            if flip:
                ntmp = n1
                n2 = n1
                n1 = ntmp
            bppm[n1, n2] += weights[i]
            if symmetric:
                bppm[n2, n1] += weights[i]
            if weight_err is not None:
                bppm_err[n1, n2] += weight_err[i]**2
                if symmetric:
                    bppm_err[n2, n1] += weight_err[i]**2

    if weight_err is not None:
        bppm_err = sqrt(bppm_err)
        for i in xrange(bppm.shape[0]):
            for j in xrange(bppm.shape[1]):
                if bppm[i, j] != 0 and bppm_err[i, j] != 0 and (bppm[i, j] / bppm_err[i, j] < signal_to_noise_cutoff):
                    bppm[i, j] = 0
        return bppm, bppm_err
    else:
        return bppm


def compare_by_register_shifts(s1, s2):
    s1bp = ss.SecondaryStructure(dbn=s1).base_pairs()
    s2bp = ss.SecondaryStructure(dbn=s2).base_pairs()
    nsharedbp = 0
    for bp1 in s1bp:
        for bp2 in s2bp:
            if (bp1[0] == bp2[0] and bp1[1] == bp2[1]) or (bp1[0] == bp2[0] - 1 and bp1[1] == bp2[1]) or (bp1[0] == bp2[0] + 1 and bp1[1] == bp2[1]) or (bp1[0] == bp2[0] + 1 and bp1[1] == bp2[1] - 1) or (bp1[0] == bp2[0] + 1 and bp1[1] == bp2[1] + 1):
                    nsharedbp += 0.5
    c = 0.8
    return nsharedbp > c * len(s1bp) and nsharedbp > c * len(s2bp)


def _collapse_by_similarity(structs):
    nstats = len(structs)
    visited = []
    state_dict = defaultdict(list)
    for i, s1 in enumerate(structs):
        if i in visited:
            continue
        for j in xrange(i+1, len(structs)):
            if j in visited:
                continue
            s2 = structs[j]
            if compare_by_register_shifts(s1, s2):
                state_dict[i].append(j)
                nstats -= 1
                visited.append(j)
    lone_structs = defaultdict(list)
    for i in xrange(len(structs)):
        found = False
        for k, v in state_dict.iteritems():
            if i in v or i == k:
                found = True
        if not found:
            lone_structs[i] = []
    state_dict.update(lone_structs)
    return state_dict


def combine_weights(w, state_dict):
    new_w = []
    for k, v in state_dict.iteritems():
        weight = w.ravel()[k]
        for s in v:
            weight += w.ravel()[s]
        new_w.append(weight)
    return array([new_w])


def classify_by_weights(weights):
    cat = 'III'
    for weight in weights.ravel():
        if weight > 0.96:
            return 'I'
        if weight > 0.5:
            cat = 'II'
    return cat


def classify(structures, weights):
    cw = combine_weights(weights, _collapse_by_similarity(structures))
    return classify_by_weights(cw)


def wentropy(weights):
    if len(weights.ravel()) == 1:
        return 0.
    else:
        return sum([-p * log(p) for p in weights.ravel() if p > 0])


def calculate_energies_from_weights(weights):
    pass

