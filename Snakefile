wildcard_constraints:
    sample="\w+",

SAMPLES, = glob_wildcards("genome_fastas/{sample}.fas")

rule all:
    input:
        "analysis/smat.bb.best.treefile",
        "reports/genomestats.csv"

rule prokka:
    input: 
       "genome_fastas/{sample}.fas"
    output:
       "annotation/{sample}/{sample}.faa",
       "annotation/{sample}/{sample}.txt",
       "annotation/{sample}/{sample}.tsv",
       "annotation/{sample}/{sample}.fsa",
    conda:
        "environment.yaml"
    threads: 8 
    shell:
        """
        prokka --outdir annotation/{wildcards.sample} --force --prefix {wildcards.sample} --locustag {wildcards.sample} --cpus {threads} --genus Wolbachia --species sp. --strain {wildcards.sample} {input}
        """

rule stats:
    input: 
        txt = expand("annotation/{sample}/{sample}.txt", sample = SAMPLES),
        tsv = expand("annotation/{sample}/{sample}.tsv", sample = SAMPLES),
        fsa = expand("annotation/{sample}/{sample}.fsa", sample = SAMPLES)
    output:
        "reports/genomestats.csv"
    conda:
        "environment.yaml"
    threads: 64 
    shell:
        """
        mkdir -p reports
        for i in {input.txt} 
        do grep 'organism\|contigs\|bases\|CDS\|rRNA\|tRNA' $i | cut -f2-20 -d' ' | sed 's/Wolbachia sp. //g' > reports/$(basename $i) 
        done
        for i in {input.fsa}
        do
        scripts/gc.sh $i >> reports/$(basename $i .fsa).txt
        scripts/cd.sh annotation/$(basename $i .fsa)/$(basename $i .fsa).tsv $i >> reports/$(basename $i .fsa).txt 
        done        
        set +e
        for j in {input.tsv}
        do 
        grep -c 'IS' $j  >> reports/$(basename $j .tsv).txt 
        grep 'CDS' $j | awk "{{sum += \$3}} END {{print sum/NR}}" | cut -f1 -d'.' >> reports/$(basename $j .tsv).txt
        done 
        set -e
        for k in reports/*.txt;
        do tr '\n' ',' < $k | tr -d ' '| sed 's/,$/\\n/' > ${{k}}1
        done
        echo "strain,contigs,size,CDS,rRNAs,tRNAs,gc,cod.dens,IS,CDS.len" > {output} 
        cat reports/*.txt1 >> {output}
        rm reports/*.txt reports/*.txt1
        checkm taxonomy_wf --tab_table -f reports/checkm.summary -t {threads} -x fas order Rickettsiales genome_fastas reports/checkm_out
        Rscript scripts/tablegen.R
        """

checkpoint orthofinder:
    input:
        expand("annotation/{sample}/{sample}.faa", sample = SAMPLES)
    output:
        OF = temp(directory("orthofinder_out")),
        SC = directory("single_copy")
    conda:
        "environment.yaml"
    threads: 64
    shell:
        """
        mkdir -p proteins
        cp {input} proteins/
        orthofinder -t {threads} -M msa -os -T iqtree -f proteins -o {output.OF}
        mkdir -p single_copy
	mv orthofinder_out/*/Single_Copy_Orthologue_Sequences/* {output.SC} 
        """

rule align:
    input:
        "single_copy/{i}.fa"
    output:
        "aligned/{i}.fas"
    conda:
        "environment.yaml"
    threads: 8
    shell:
        """     
        mkdir -p aligned 
        linsi --thread {threads} {input} | sed 's/_.*//g' > {output}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.orthofinder.get(**wildcards).output[1]
    return expand('aligned/{i}.fas',
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.fa')).i)

rule concat:
    input: 
        aggregate_input
    output:
        matrix= "analysis/supermat.fas",
        partition= "analysis/partition.txt" 
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p analysis
        mkdir -p aligned/recombining
        for j in {{10,20,30,40,50,100}}
        do
          for i in {input}
          do 
          Phi  -f ${{i}} -w ${{j}} -t A > $(basename ${{i}} .fas).philog
          done
        grep 'PHI (Normal):' *.philog | cut -f1,3 -d':' | sed 's/.philog://' > aligned/recombining/recombination_w${{j}}.res
        sort -g -k2 aligned/recombining/recombination_w${{j}}.res | awk '$2 < 0.05 && $2!= "--"' | cut -f1 -d' ' | xargs -I{{}} cp aligned/{{}}.fas aligned/recombining/
        rm *.philog Phi*
        done
	ls aligned/recombining/ | grep '.fas' | xargs -I{{}} rm aligned/{{}}
        pxcat -p {output.partition} -o {output.matrix} -s aligned/*.fas
        rm phyx.logfile
        """

rule iqtree:
    input: 
        "analysis/supermat.fas"
    output:
        "analysis/smat.fast.LG.treefile",
        "analysis/smat.bb.best.treefile"
    conda:
        "environment.yaml"
    threads: 48
    shell:
        """
        iqtree -m LG+I+G4 -s {input} -pre analysis/smat.fast.LG -redo -fast -nt {threads}
        iqtree -s {input} -pre analysis/smat.bb.best -bb 1000 -redo -nt {threads}
        """

