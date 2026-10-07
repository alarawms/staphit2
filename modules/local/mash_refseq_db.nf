process MASH_REFSEQ_DB {
    label 'process_low'
    container 'ghcr.io/alarawms/staphit2-python:1.0.0'

    output:
    path "refseq.genomes.msh", emit: db
    tuple val("${task.process}"), val('mash_refseq_md5'), eval("md5sum refseq.genomes.msh | cut -d' ' -f1"), topic: versions, emit: versions_mash_refseq_md5

    script:
    '''
    python3 -c "
import urllib.request
urllib.request.urlretrieve(
    'https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh',
    'refseq.genomes.msh'
)
print('Downloaded RefSeq Mash sketch')
"
    '''

    stub:
    """
    touch refseq.genomes.msh
    """
}
