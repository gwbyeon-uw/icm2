import os
# import pdb

import matplotlib
from matplotlib.pylab import *
from matplotlib.colors import rgb2hex
from matplotlib.patches import Rectangle

from rdatkit import VARNA, SecondaryStructure

import map_analysis_utils as utils

STRUCTURE_COLORS = [get_cmap('Paired')(i * 50) for i in xrange(100)]


def plot_mutxpos_image(d, sequence, seqpos, offset, mut_labels, cmap=get_cmap('Greys'), vmin=0, vmax=None, missed_indices=None, contact_sites=None, structure_colors=STRUCTURE_COLORS, weights=None, aspect='auto'):
    ax = subplot(111)
    if vmax is None:
        vmax = d.mean()
    ax.imshow(d, cmap=get_cmap('Greys'), vmin=0, vmax=vmax, aspect=aspect, interpolation='nearest')
    if missed_indices is not None:
        for x, y in missed_indices:
            ax.add_artist(Rectangle(xy=(y - 0.5, x - 0.5), facecolor='none', edgecolor='r', linewidth=0.25, width=1, height=1))
    if contact_sites is not None:
        if weights is None:
            weights = ones(d.shape[0], len(contact_sites))
        for k, v in contact_sites.iteritems():
            for x, y in zip(*where(v != 0)):
                ax.add_artist(Rectangle(xy=(y - 0.5, x - 0.5), facecolor='none', edgecolor=structure_colors[k], linewidth=1, width=1, height=1, alpha=weights[x, k]))
    xticks(range(len(seqpos)), ['%s%s' % (pos, sequence[pos - offset - 1]) for pos in seqpos], fontsize='xx-small', rotation=90)
    yticks(range(len(mut_labels)), mut_labels, fontsize='xx-small')
    return ax


def expected_reactivity_plot(react, struct, yerr=None, ymin=0, ymax=5, seq_indices=None):
    if seq_indices is None:
        seq_indices = [i for i in xrange(len(react))]
    binarized_react = [1 if x == '.' else 0 for i, x in enumerate(struct[:len(react)]) if i in seq_indices]
    plot(1 + arange(len(binarized_react)), binarized_react, linewidth=2, c='r')
    bar(0.5 + arange(len(react)), react, linewidth=0, width=1, color='gray')
    if yerr is not None:
        errorbar(1 + arange(len(react)), react, yerr=yerr, linewidth=2, c='k')
    ylim(0, 5)
    xlim(1, len(react))


def weights_by_mutant_plot(W, W_err, mut_labels, structure_colors=STRUCTURE_COLORS, W_ref=None, W_samples=None, idx=-1, assignments=None, medoids=None):
    ax = subplot(111)
    if assignments is None:
        _W = W
        _W_err = W_err
        _W_ref = W_ref
        nstructs = W.shape[1]
        struct_indices = range(nstructs)
    else:
        _W = zeros([W.shape[0], len(medoids)])
        _W_err = zeros([W.shape[0], len(medoids)])
        if W_ref is not None:
            _W_ref = zeros([W.shape[0], len(medoids)])
        else:
            _W_ref = None
        nstructs = len(medoids)
        struct_indices = medoids
        i = 0
        setidx = False
        for c, si in assignments.iteritems():
            if idx in si:
                if not setidx:
                    idx = i
                    setidx = True

            _W[:, i] = W[:, si].sum(axis=1)

            if _W_ref is not None:
                _W_ref[:, i] = W_ref[:, si].sum(axis=1)

            if W_samples is not None:
                """
                a = ones(len(si))
                for j in xrange(_W_err.shape[0]):
                    _W_err[j,i] = dot(dot(a, cov(W_samples[j,si,:])), a.T)
                print sqrt(_W_err[0,:])
                """
                _W_err[:, i] = W_samples[:, si, :].sum(axis=1).std(axis=1)
            else:
                _W_err[:, i] = (W_err[:, si]**2).sum(axis=1)
            #_W_err[:,i] = sqrt(_W_err[:,i])

            i += 1
    if idx >= 0:
        weight_range = [idx]
    else:
        weight_range = xrange(nstructs)

    for j in weight_range:
        ax.errorbar(arange(W.shape[0])+1, _W[:, j], yerr=_W_err[:, j], linewidth=3, label='structure %s ' % (struct_indices[j]), color=structure_colors[j])
        if _W_ref is not None:
            ax.errorbar(arange(_W_ref.shape[0]) + 1, _W_ref[:, j], linestyle='--', linewidth=3, label='reference %s ' % (struct_indices[j]), color=structure_colors[j])

    ylim(0, 1)
    xlim(0, _W.shape[0] + 1)
    xticks(arange(_W.shape[0]) + 1, mut_labels, fontsize='xx-small', rotation=90)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return _W, _W_err, _W_ref


def PCA_structure_plot(structures, assignments, medoids, colorbyweight=False, weights=None, names=None):
    all_struct_vecs, all_struct_indices, select_struct_vecs = [], [], []
    cluster_colors, struct_to_clust = {}, {}

    for i, c in enumerate(assignments):
        cluster_colors[c] = STRUCTURE_COLORS[i]
    for c, indices in assignments.iteritems():
        for i in indices:
            struct_to_clust[i] = c
            all_struct_indices.append(i)
            vec = [1 if s == '.' else 0 for s in structures[i]]
            all_struct_vecs.append(vec)
            if i in medoids:
                select_struct_vecs.append(vec)
    all_struct_vecs = array(all_struct_vecs)
    select_struct_vecs = array(select_struct_vecs)

    U, s, Vh = svd(all_struct_vecs.T)
    basis = U[:, :2]
    all_struct_coordinates = dot(all_struct_vecs, basis)
    select_struct_coordinates = dot(select_struct_vecs, basis)

    all_sizes, medoid_sizes = 50, 100
    if weights is None:
        all_structure_colors = [cluster_colors[struct_to_clust[i]] for i in all_struct_indices]
        medoid_colors = [cluster_colors[struct_to_clust[i]] for i in medoids]
        linewidth = 1
    else:
        cmap = get_cmap('jet')
        normf = Normalize(vmin=0, vmax=max(weights), clip=True)
        if colorbyweight:
            all_structure_colors = [cmap(normf(weights[i])) for i in all_struct_indices]
            medoid_colors = [cmap(normf(weights[i])) for i in medoids]
            linewidth = 0
        else:
            all_structure_colors = [cluster_colors[struct_to_clust[i]] for i in all_struct_indices]
            medoid_colors = [cluster_colors[struct_to_clust[i]] for i in medoids]
            all_sizes = array([max(50, log(1 + weights[i])*5e4) for i in all_struct_indices])
            medoid_sizes = array([max(50, log(1 + weights[i])*5e4) for i in medoids])
            linewidth = 1

    figure(1)
    clf()
    scatter(all_struct_coordinates[:, 0], all_struct_coordinates[:, 1], c=all_structure_colors, alpha=0.6, linewidth=linewidth, s=all_sizes)
    scatter(select_struct_coordinates[:, 0], select_struct_coordinates[:, 1], c=medoid_colors, linewidth=2, s=medoid_sizes)
    if names is not None:
        for i, m in enumerate(medoids):
            if '_' in names[i]:
                names[i] = names[i].replace('_', '_{') + '}'
            text(select_struct_coordinates[i, 0], select_struct_coordinates[i, 1], '$' + names[i] + '$', style='italic')


def bpp_matrix_plot(structures, weights, ref_weights=None, weight_err=None, offset=0):
    if weight_err is not None:
        bppm, bppm_err = utils.bpp_matrix_from_structures(structures, weights, weight_err=weight_err, signal_to_noise_cutoff=1)
    else:
        bppm = utils.bpp_matrix_from_structures(structures, weights)

    if ref_weights is not None:
        bppm_ref = utils.bpp_matrix_from_structures(structures, ref_weights)
        for i in xrange(bppm_ref.shape[0]):
            for j in xrange(i + 1, bppm_ref.shape[1]):
                bppm[i, j] = bppm_ref[j, i]

    r = arange(0, bppm.shape[0], 10)
    r_offset = int(offset / 10) * 10 + 10 - offset - 1
    r = array([0] + (r + r_offset).tolist())
    r[-1] = min(bppm.shape[0] - 1, r[-1])

    colors = [(cm.jet(i)) for i in xrange(235)]
    bppm_map = matplotlib.colors.LinearSegmentedColormap.from_list('bppm_map', colors)
    #bppm[bppm <= 0.05] = 0
    imshow(bppm, cmap=bppm_map, interpolation='nearest', vmax=1)
    grid()
    colorbar()
    xticks(r, r + offset + 1, rotation=90)
    yticks(r, r + offset + 1)
    return bppm


def make_struct_figs(structures, sequence, offset, fprefix, indices=None, base_annotations=None, helix_function=lambda x,y:x, helix_fractions=None, annotation_color='#FF0000'):
    options = {'drawBases': False, 'fillBases': False, 'resolution': '10.0', 'flat': True, 'offset': offset}
    if indices is None:
        indices = range(len(structures))
    for i, s in enumerate(structures):
        print s
        options['bp'] = rgb2hex(STRUCTURE_COLORS[i])
        varna = VARNA(sequences=[sequence], structures=[SecondaryStructure(dbn=s)])
        if base_annotations is None:
            CMD = varna.render(output=fprefix + 'structure%s.svg' % indices[i], annotation_by_helix=True, helix_function=helix_function, cmd_options=options)
        else:
            varna.annotation_font_size = 13
            varna.annotation_color = annotation_color
            if helix_fractions is None:
                helix_frac_annotations = ''
                base_weight_annotations = varna._get_base_annotation_string([base_annotations[i]], annotation_by_helix=True, helix_function=helix_function)
            else:
                helix_frac_annotations = varna._get_base_annotation_string([helix_fractions[i]], annotation_by_helix=True, helix_function=helix_function, stype='B', helix_side=0)
                #varna.annotation_color = '#0033CC'
                #base_weight_annotations = varna._get_base_annotation_string([base_annotations[i]], annotation_by_helix=True, helix_function=helix_function, stype='B', helix_side=0, base_offset=-2)
            options['annotations'] = helix_frac_annotations.strip('"')
            #options['annotations'] += base_weight_annotations.strip('"')
            CMD = varna.render(output=fprefix + 'structure%s.svg' % indices[i], annotation_by_helix=True, helix_function=helix_function, cmd_options=options)
        print CMD
        os.system(CMD)

