process MASH_REFSEQ_DB {
    label 'process_low'
    publishDir "${params.outdir}/databases", mode: 'copy'
    container 'docker.io/python:3.9'

    output:
    path "refseq.genomes.msh", emit: db

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
}
