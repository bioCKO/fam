#!/usr/bin/env python3

import os
import sys
import seq
import subprocess
import re

syn_id_file = sys.argv[1]
group_file = sys.argv[2]
cds_file= sys.argv[3]
pep_file= sys.argv[4]
out_dir = sys.argv[5]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
else:
    print(out_dir + ' exists!')
    sys.exit()

syn_genes = set()
for line in open(syn_id_file, 'rt'):
    line = line.strip()
    if line:
        ss = line.split('\t')
        if len(ss) >= 2:
            syn_genes.add(ss[0])
            syn_genes.add(ss[1])

group2genes = {}
gene2group = {}
group2syn_genes = {}
for line in open(group_file, 'rt'):
    line = line.strip()
    if line:
        ss = line.split('\t')
        if len(ss) == 2:
            group_id = ss[0]
            gene_id = ss[1]
            if group_id in group2genes:
                group2genes[group_id].append(gene_id)
            else:
                group2genes[group_id] = [gene_id]
            gene2group[gene_id] = group_id

            if gene_id in syn_genes:
                if group_id in group2syn_genes:
                    group2syn_genes[group_id].append(gene_id)
                else:
                    group2syn_genes[group_id] = [gene_id]

gene_names, name2cds = seq.fasta_seq(cds_file)
gene_names, name2pep = seq.fasta_seq(pep_file)

out_file = os.path.join(out_dir, 'out.txt')
out = open(out_file, 'wt')

for group_id in group2syn_genes:
    gene_ids = group2genes[group_id]
    syn_genes = group2syn_genes[group_id]
    group_dir = os.path.join(out_dir, group_id)
   
    cds_file = os.path.join(group_dir, 'cds.fa')        
    cds_aln_file = os.path.join(group_dir,'cds_aln.fa')
    no_gap_cds_aln_file = os.path.join(group_dir, 'no_gap_cds_aln.fa')
    paml_aln_file = os.path.join(group_dir, 'paml.aln')
    pep_file = os.path.join(group_dir, 'pep.fa')
    pep_aln_file = os.path.join(group_dir,'pep_aln.fa')
    no_gap_pep_aln_file = os.path.join(group_dir, 'no_gap_pep_aln.fa')
    tree_file = os.path.join(group_dir, 'ml.nwk')
    
    if not os.path.exists(group_dir):
        os.mkdir(group_dir)
    
    fout = open(cds_file, 'wt')
    for gene_id in gene_ids:
        if gene_id in name2cds:
            fout.write('>' + gene_id + '\n')
            fout.write(name2cds[gene_id] + '\n')
        else:
            print("No cds for: " + gene_id)
            sys.exit()
    fout.close()

    seq.aln_cds(cds_file, cds_aln_file)
    seq.remove_aln_gap(0, cds_aln_file, no_gap_cds_aln_file)

    gene_names, gene2seq= seq.fasta_seq(no_gap_cds_aln_file)
    fout = open(paml_aln_file, 'wt')
    s = ' {} {}\n\n'.format(len(gene_names), len(gene2seq[gene_names[0]]))
    fout.write(s)
    for gene_name in gene_names:
        gene_seq = gene2seq[gene_name]
        fout.write(gene_name + '\n')
        fout.write(gene_seq + '\n\n')
    fout.close()

    fout = open(pep_file, 'wt')
    for gene_id in gene_ids:
        if gene_id in name2pep:
            fout.write('>' + gene_id + '\n')
            fout.write(name2pep[gene_id] + '\n')
        else:
            print("No pep for: " + gene_id)
            sys.exit()
    fout.close()

    com = 'muscle -quiet -in {} -out {}'.format(pep_file, pep_aln_file)
    subprocess.call(com, shell=True)
    seq.remove_aln_gap(0, pep_aln_file, no_gap_pep_aln_file)
    com = 'fasttree -wag -gamma -nosupport -out {} {} '.format(tree_file, no_gap_pep_aln_file)
    subprocess.call(com, shell=True)

    for syn_gene in syn_genes:
        if len(gene_ids) < 3:
            continue
        paml_tree_file = os.path.join(group_dir, syn_gene + '.nwk')
        mlc_null_file = os.path.join(group_dir, syn_gene + '_mlc_null')
        mlc_alter_file = os.path.join(group_dir,syn_gene +  '_mlc_alter')
        ctl_null_file = os.path.join(group_dir, syn_gene + '_null.ctl')
        ctl_alter_file = os.path.join(group_dir, syn_gene + '_alter.ctl')
        p_file = os.path.join(group_dir, syn_gene + '_p.txt')

        nwk = ''
        for line in open(tree_file, 'rt'):
            line = line.strip()
            if line:
                nwk += line

        nwk = nwk.replace(syn_gene, syn_gene+' #1 ')
        fout = open(paml_tree_file, 'wt')
        fout.write(str(len(gene_ids)) + ' 1\n\n')
        fout.write(nwk + '\n')
        fout.close()

        fout = open(ctl_null_file, 'wt')
        s = '''
        seqfile = {}
     treefile = {}
      outfile = {}
        noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
   aaRatefile = wag.dat
        model = 2
      NSsites = 2
        icode = 0
    fix_kappa = 0
        kappa = 3
    fix_omega = 1
        omega = 1
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 10
        getSE = 0
 RateAncestor = 0
    Small_Diff = .5e-6
        '''.format(paml_aln_file, paml_tree_file, mlc_null_file)
        fout.write(s)
        fout.close()

        fout = open(ctl_alter_file, 'wt')
        s = '''
        seqfile = {}
     treefile = {}
      outfile = {}
        noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
   aaRatefile = wag.dat
        model = 2
      NSsites = 2
        icode = 0
    fix_kappa = 0
        kappa = 3
    fix_omega = 0
        omega = 1.5
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 10
        getSE = 0
 RateAncestor = 0
    Small_Diff = .5e-6
        '''.format(paml_aln_file, paml_tree_file, mlc_alter_file)
        fout.write(s)
        fout.close()

        com = 'codeml ' + ctl_null_file
        subprocess.call(com, shell=True)
        com = 'codeml ' + ctl_alter_file
        subprocess.call(com, shell=True)

        v_null = 0
        v_alter = 0
        for line in open(mlc_null_file, 'rt'):
            line = line.strip()
            if line:
                m = re.search(r'lnL.+\):\s+([\S]+)', line)
                if m:
                    v_null = float(m.group(1))  
        for line in open(mlc_alter_file, 'rt'):
            line = line.strip()
            if line:
                m = re.search(r'lnL.+\):\s+([\S]+)', line)
                if m:
                    v_alter = float(m.group(1))
        
        detal2 = 2 * abs(v_alter - v_null)
        os.system('chi2 1 {} > {}'.format(detal2, p_file))
        p_v = None
        for line in open(p_file, 'rt'):
            line = line.strip()
            if line:
                m = re.search(r'prob = ([\d\.]+) =', line)
                if m:
                    p_v = m.group(1)
        s = syn_gene + '\t' + p_v + '\n'
        print(s)
        out.write(s)
out.close()
   