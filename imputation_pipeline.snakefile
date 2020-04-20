fullcov = [20.009, 22.0742, 22.2528, 22.8935, 22.3157, 22.1753, 21.8395, 21.905, 18.9248, 21.1734, 21.3123, 21.6016, 19.0814, 18.3776, 17.4509, 18.6179, 19.5518, 21.57, 18.5276, 19.802, 17.3075, 13.4971]

rule target:
    input:
        expand('{cov}x/impute_out/NE1.chr{chrom}.atlas.{cov}X_MaximumLikelihood.gl.gp0.99.imputed.gp0.99.vcf.gz', chrom = range(1, 23), cov = [0.05, 0.1, 0.5, 0.75, 1, 1.5, 2])

rule downsample:
    input:
        inbam = 'NE1/NE1.{chrom}.bam'
    output:
        outbam = 'downsampled_bams/NE1.chr{chrom}.cov{cov}X.bam'
    params:
        frac = lambda wildcards: float(wildcards.cov) / fullcov[int(wildcards.chrom) - 1] + 10
    shell:
        '''
        samtools view -s {params.frac} -b {input.inbam} > {output.outbam}
        samtools index {output.outbam}
        '''

rule atlas_call:
    input: 'downsampled_bams/NE1.chr{chrom}.cov{cov}X.bam'
    output: 'atlas/NE1.chr{chrom}.atlas.{cov}X_MaximumLikelihood.vcf.gz'
    params:
        outname = lambda wildcards, output: output[0].replace('_MaximumLikelihood.vcf.gz', '')
    shell:
        '''
        sites=1kg_sites/1kg.chr{wildcards.chrom}.sites
        atlas task=call method=MLE infoFields=DP formatFields=GT,AD,DP,PL bam={input} chr={wildcards.chrom} fasta=resources/hs37d5/hs37d5.fa sites=${{sites}} out={params.outname}
        '''

rule gp_filter:
    input:
        infile='{filestem}.vcf.gz'
    output:
        outfile='{filestem}.gp0.99.vcf.gz'
    shell:
        '''
        bcftools view -i 'MAX(GP)>=0.99' {input.infile} -Oz -o {output.outfile}
        '''

rule beagle_gl:
    input:
        infile='atlas/NE1.chr{chrom}.atlas.{cov}X_MaximumLikelihood.vcf.gz',
        glreffile='imputation_reference_files/chr{chrom}.1kg.phase3.v5a.vcf.gz',
        mapfile='genetic_mapGRCh37/plink.chr{chrom}.GRCh37.map'
    output:
        '{cov}x/gl_out/NE1.chr{chrom}.atlas.{cov}X_MaximumLikelihood.gl.vcf.gz'
    params:
        outname = lambda wildcards, output: output[0].replace('.vcf.gz', '')
    shell:
        '''
	java -Xss10m -Xmx16g -jar beagle.11Mar19.69c.jar gl={input.infile} ref={input.glreffile} map={input.mapfile} out={params.outname} nthreads=${{SLURM_CPUS_PER_TASK}} gprobs=true chrom={wildcards.chrom}
        '''
        
rule beagle_gt:
    input:
        infile='{cov}x/gl_out/NE1.chr{chrom}.atlas.{cov}X_MaximumLikelihood.gl.gp0.99.vcf.gz',
        mapfile='genetic_mapGRCh37/plink.chr{chrom}.GRCh37.map',
        gtreffile='imputation_reference_files/chr{chrom}.1kg.phase3.v5a.vcf.gz'
    output:
        '{cov}x/impute_out/NE1.chr{chrom}.atlas.{cov}X_MaximumLikelihood.gl.gp0.99.imputed.vcf.gz'
    params:
        outname = lambda wildcards, output: output[0].replace('.vcf.gz', '')
    shell:
        '''
        java -Xss10m -Xmx32g -jar beagle.16May19.351.jar gt={input.infile} ref={input.gtreffile} out={params.outname} impute=true nthreads=${{SLURM_CPUS_PER_TASK}} map={input.mapfile} gp=true
        '''
