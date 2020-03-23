###########################################################################
## Make figures and compile methods tex file for paper 
#########################################################################

shell.prefix("source ~/.bashrc; ")


configfile: "config.yaml"

import pandas as pd
import subprocess
import os
import math              
            
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()            
            
localrules: all

import os

subworkflow GeuRefbias:
    workdir: "/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias"

subworkflow GeuRefbias2:
    workdir: "/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias2"             

subworkflow InputPrep:
    workdir: "/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis"

subworkflow PsoBtrecase:
    workdir: "/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase"

subworkflow PsoRefBias:
    workdir: "/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/refbias"
             


# get home directories for Snakefiles in subdirectories to call Scripts/Functions etc avoiding duplication
home_GeuRefbias = vars(vars(workflow)['_subworkflows']['GeuRefbias'])['_workdir']
home_GeuRefbias2 = vars(vars(workflow)['_subworkflows']['GeuRefbias2'])['_workdir']
home_InputPrep = vars(vars(workflow)['_subworkflows']['InputPrep'])['_workdir']
home_PsoBtrecase = vars(vars(workflow)['_subworkflows']['PsoBtrecase'])['_workdir']
home_PsoRefBias = vars(vars(workflow)['_subworkflows']['PsoRefBias'])['_workdir']


## These functions are from /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Snakefile

def gene_chrom(File=PsoBtrecase(config['pso_dir']) + "/Btrecase/inputs/gene_inputs/genes2test.txt", sep=" "):
    """ Makes a dictionary with keys gene_id and values chromosome from a file with first col gene_id and second col CHROM """
    data=pd.read_csv(File, sep=sep)
    keys=list(data['gene_id'])
    values=[str(x) for x in data['CHROM']]
    dic=dict(zip(keys,values))
    return dic

def geneName_geneID(File=config['pso_dir'] + "/Btrecase/results/full.stan.summary.txt", sep=" "):
    """ Makes a dictionary with keys gene_name and values gene_id """
    data=pd.read_csv(File, sep=sep)
    keys=list(data['Gene_name'])
    values=[str(x) for x in data['Gene_id']]
    dic=dict(zip(keys,values))
    return dic

pso_genes=["ERAP1", "PPIF"]
pso_id=[geneName_geneID()[g] for g in pso_genes]  ## gene_id
chrom_pso=[gene_chrom()[c] for c in pso_id]  ## chrom

rule all:
    input:
        #expand(config['out_dir'] + "/{IndVar}.png", IndVar=["Between", "Within"])
        #expand(config['out_dir'] + "/Figures/{fig}.pdf", fig=["Fig" + str(x) for x in range(1,2)])
        #expand("{fig}.pdf", fig=["Fig" + str(x) for x in range(1,2)] )
        #expand(config['out_dir'] + "/Figures/{fig}.png", fig=["Fig5"])
        expand(config['var_dir'] + "/variation_homo_sapiens-chr{chrom}.gvf.gz", chrom=[22] + chrom_pso)


rule graphs_fig1:
    """ Prepare graphs to include in Figure1 """
    output:
        outf=expand(config['out_dir'] + "/Figures/{IndVar}.png", IndVar=["Between", "Within"])
    script:
        "Scripts/between.ind.gene.R"

rule down_SNPs:
    """Download file with variant information by chromosomes, alternative to biomart"""
    input:
        FTP.remote("ftp.ensembl.org/pub/grch37/current/variation/gvf/homo_sapiens/homo_sapiens-chr{chrom}.gvf.gz",
                   keep_local=True, immediate_close=True)
    output:
        config['var_dir'] + "/variation_homo_sapiens-chr{chrom}.gvf.gz"
    run:
        outputName =os.path.join(config['var_dir'] ,"variation_" + os.path.basename(input[0]))
        shell("mv {input} {outputName}" )
 
        
rule geu_figs:
    """ Prepares figures from geuvadis data"""
    input:
        geneStEnd=InputPrep(config['geneStEnd']), 
        AIgt=GeuRefbias2(config["geu_refb2"] + "/post_remap/pre_post_AI.txt"),
        AIngt=GeuRefbias(config["geu_refb"] + "/post_remap/pre_post_AI.txt"),
        eSNPs=config['fSNPs'],
        dseq=expand(config['dseq'] + "/RunNBmodelbatch{N}_chr22.nbmodelRes.csv", N=[x+1 for x in range(10)]),
        gtex= config['EBV'] + "/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz",
        sigGtex=config['EBV'] + "/Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz",
        #geu_eur=config['EUR-GEU'],
        geu_chris=config["sig-GEU"],
        geu_egenes=GeuRefbias(config["geu_arrayX"]) + "/EUR373.gene.cis.FDR5.best.rs137.txt.gz"
        rna_dp=config["rna_qc"],
        var22=config['var_dir'] + "/variation_homo_sapiens-chr22.gvf.gz",
        het_prop=GeuRefbias(config['geu_refb']) + "/Btrecase/QCfSNP/fisher.fSNPs.txt",
        dna=config['dna'],
        rna=config['rna'],
        #rna_dp=config['rna_depth'],
        sample_down_date=config['geu'] + "sample_info/sample.downdate",
        sample_info=config['geu'] + "sample_info/GBR.sample"
    params:
        #btrecase=expand(GeuRefbias(config['geu_refb']) + "/Btrecase/{prior}/{source}", prior=["SpikeMixV3_2"], source=["GT","RNA"]),
        btrecase=expand(GeuRefbias(config['geu_refb']) + "/Btrecase/SpikeMixV3_2/{source}", source=["GT", "fisher001/RNA"]),
        btrecaseGC=GeuRefbias(config['geu_refb2']) + "/Btrecase/SpikeMixV3_2/GT"
        lm=GeuRefbias(config['geu_refb']) + "/lm/log_counts",
        trec=expand(GeuRefbias(config['geu_refb']) + "/Btrec/{prior}", prior=["SpikeMixV3_2"]),
        rbias=["norefbias", "rbias"],
        rasqual=GeuRefbias(config['geu_refb']) + "/rasqual/outCounts/cis5_10_5",
        rasqual01=GeuRefbias(config['geu_refb']) + "/rasqual/output/cis1_10_5",  ## rasqual smaller cis-window
        QC2=GeuRefbias(config['geu_refb']) + "/Btrecase/QC2_F/RNA",
        QC1=GeuRefbias(config['geu_refb']) + "/Btrecase/QC1_F/RNA",
        flatPrior=GeuRefbias(config['geu_refb']) + "/Btrecase/weakprior/GT",
        sim_dir=config['simdir'],
        sim_pattern=["diff_0.05_hfq_0.1_" + x for x in ["pop", "beta"]], ## simulation files
        inf_prior=[0,0, 0.0309, 0.3479, 0.97359164, 0.02640836], # mean, sd and mixing proportion,
        non_infprior=[0,10, 1], #mean and sd and mix
        out=config['out_dir'] + "/Figures"
    output:
        expand(config['out_dir'] + "/Figures/{fig}.pdf", fig=["Fig2bottom", "Fig3"])
    script:
        "Scripts/geu_figs.R"


rule latex_figs:
    """ Prepare latex figures, diagram with method description. Coudn't get output-dir command in pdflatex to work from snakemake. Also, I added auxfig in input to link rules but the graphs are called manually in latex. Also graphicspath was set manually in document """
    input:
         Fig="Scripts/{fig}.tex",
         auxfig=expand(config['out_dir'] + "/Figures/{IndVar}.png", IndVar=["Between", "Within", "Fig2bottom"])
    params:
        out=config['out_dir']
    output:
        config['out_dir'] + "/Figures/{fig}.pdf"
        #"{fig}.pdf"
    shell:
        "pdflatex  {input.Fig} ; "
        "rm *aux ;"
        "rm *log ;"
        "mv {wildcards.fig}.pdf {output} "
  
        
rule nor_pso:
    """ Prepare figures running normal and psoriasis in different models """
    input:
        #oldsum=PsoBtrecase(config['pso_dir']) + "/Btrecase/results/normal_pso_ci99_summary.txt",
        gene_coord=InputPrep(config['pso_dir']) + "/Btrecase/inputs/gene_inputs/gene_info.txt",
        gtexSigExp=config['EBV'] + "/Skin_Sun_Exposed_Lower_leg.v7.signif_variant_gene_pairs.txt.gz",
        gtexSigNoexp=config['EBV'] + "/Skin_Not_Sun_Exposed_Suprapubic.v7.signif_variant_gene_pairs.txt.gz",
        var=expand(config['var_dir'] + "/variation_homo_sapiens-chr{chrom}.gvf.gz", chrom=chrom_pso),
        ##r2=PsoRefBias(config['pso_dir']) + "/Btrecase/objects/EURr2.rds",
        drg=PsoBtrecase(config['pso_dir']) + "/pso_gwas/DRG_JID2014.xls",
        gwas=PsoBtrecase(config['pso_dir']) + "/pso_gwas/SupTable2.csv",
        cor=PsoRefBias(config['pso_dir']) + "/Btrecase/objects/fisher001EURr.rds",
        gwasr2=PsoRefBias(config['pso_refbias_dir']) + "/Btrecase/results/r2_spike_gwas.txt",
    params:
        stan_dir=PsoRefBias(config['pso_refbias_dir']) + "/Btrecase/SpikePrior/fisher001/rna",
        stan1M2T_dir=PsoRefBias(config['pso_refbias_dir']) + "/Btrecase/2Tissues/jointInfo_Fish001_hetInd",
        gtex_dir=config['pso_dir'] + "/eQTL_gtex",
        skin=["normal_skin", "Psoriasis_skin"],
        d=10**5,
        infoMax=1.1,
        genes2follow=["ERAP1", "PPIF"],
        rbias="refbias", # refbias correction only
        y_G1=[1.4,2.4], # genes2follow 1
        y_G2=[1.4, 1.6, 1.4, 2.4, 1.4], # genes2follow 2
        j_G1=[0]*2 , #list with just values G1
        j_G2=[0]*5,
        r_G1=[1.7, 2.3], # rect for G1 (y axis)
        r_G2=[1.7, 2.3], # rect for G2
        w=3000,
        colclass=["character"]+["numeric" for i in range(7)]
    output:
        expand(config['out_dir'] + "/Figures/{fig}.png", fig=["Fig5", "SuppTable3.csv"])
    script:
        "Scripts/nor_pso.R"
        
