process FETCH_RESFINDER_DB {
    label 'process_low'
    container 'docker.io/python:3.9'

    output:
    path "resfinder_db", emit: db

    script:
    '''
    mkdir -p resfinder_db
    python3 << 'PYEOF'
import urllib.request, zipfile, io, os
url = 'https://bitbucket.org/genomicepidemiology/resfinder_db/get/master.zip'
resp = urllib.request.urlopen(url)
z = zipfile.ZipFile(io.BytesIO(resp.read()))
count = 0
for f in z.namelist():
    if f.endswith('.fsa'):
        outname = os.path.basename(f)
        with open(os.path.join('resfinder_db', outname), 'wb') as out:
            out.write(z.read(f))
        count += 1
print(f'Downloaded {count} ResFinder database files')
PYEOF
    '''
}

process INDEX_DB {
    label 'process_low'
    container 'docker.io/staphb/kma:1.4.14'

    input:
    path db

    output:
    path "indexed_db", emit: indexed_db

    script:
    """
    mkdir indexed_db
    kma index -i ${db}/*.fsa -o indexed_db/resfinder
    """
}

process KMA {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/staphb/kma:1.4.14'

    input:
    tuple val(meta), path(reads)
    path indexed_db

    output:
    tuple val(meta), path("*.res"), emit: results

    script:
    def nano = meta.mode == 'long' ? '-bcNano' : ''   // long-only samples get raw ONT reads
    """
    kma -i ${[reads].flatten().join(' ')} -o ${meta.id} -t_db indexed_db/resfinder -1t1 ${nano} || true
    touch ${meta.id}.res
    """
}
