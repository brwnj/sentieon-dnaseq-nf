nextflow.enable.dsl=2

params.bams = false
params.fasta = false
params.dbsnp = false

fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
dbsnp = file(params.dbsnp)
dbsnpidx = file("${params.dbsnp}.tbi")

bams_ch = Channel.fromPath(params.bams, checkIfExists: true)
    // sample, bam, bai, recal
    .map { file -> tuple(file.simpleName, file, file + '.bai', file.toString().replace('.deduped.bam', '.recal.table')) }

process dnaseq {
    publishDir path: "results/dnaseq"
    cpus 8

    input:
    tuple(val(sample), path(bam), path(bai), path(recal))
    path(fasta)
    path(faidx)
    path(dbsnp)
    path(dbsnpidx)

    output:
    path("${sample}_dnaseq.g.vcf"), emit: vcf

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${bam} \
        -q ${recal} \
        --algo Haplotyper \
        --emit_mode gvcf \
        -d ${dbsnp} \
        ${sample}_dnaseq.g.vcf
    """
}

process vcfidx {
    publishDir path: "results/dnaseq"

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.g.vcf.idx"), emit: idx

    script:
    """
    sentieon util vcfindex ${vcf}
    """
}

process gvcftyper {
    publishDir path: "results/gvcftyper"
    cpus 16

    input:
    path(gvcfs)
    path(gvcfidxs)
    path(fasta)
    path(faidx)
    path(dbsnp)
    path(dbsnpidx)

    output:
    path("output-joint.vcf.gz")
    path("output-joint.vcf.gz.tbi")

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        --algo GVCFtyper \
        output-joint.vcf \
        ${gvcfs}
    bgzip output-joint.vcf
    tabix -p vcf output-joint.vcf.gz
    """
}

workflow {
    dnaseq(bams_ch, fasta, faidx, dbsnp, dbsnpidx)
    vcfidx(dnaseq.output.vcf)
    gvcftyper(dnaseq.output.vcf.collect(), vcfidx.output.idx.collect(), fasta, faidx, dbsnp, dbsnpidx)
}
