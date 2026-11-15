# SPDX-License-Identifier: MIT
"""LosslessSegmentCASV1: byte-preserving, owned-segment payload retention.

This is an external output controller, not a solver cache. No scientific check,
state or history is removed. Only newly produced verified payloads are pruned;
native reports and native restart histories stay directly readable.
"""
import copy,gzip,hashlib,json,os,shutil,time
from pathlib import Path
from rans_request_manifest import sha
from ct_accepted_feasibility import require,write_json


def decoded(blob,destination=None):
    h=hashlib.sha256();size=0
    stream=gzip.open(blob,'rb');out=None
    try:
        if destination is not None:
            destination.parent.mkdir(parents=True,exist_ok=True);out=destination.open('xb')
        for block in iter(lambda:stream.read(1024*1024),b''):
            h.update(block);size+=len(block)
            if out is not None:out.write(block)
    finally:
        stream.close()
        if out is not None:out.close()
    return h.hexdigest(),size


def validate(root,archive,index,proof):
    require(index['schema']=='LosslessSegmentCASV1' and proof['status']=='pass','unqualified retention request')
    require(proof['binding']==index['binding'],'foreign retention/native identity')
    for name,digest in proof['required_files'].items():
        file=root/name;require(file.is_relative_to(root) and sha(file)==digest,'required milestone/history evidence changed:'+name)
    for name,row in index['files'].items():
        file=root/name;require(file.is_relative_to(root) and '..' not in Path(name).parts,'foreign prune path')
        blob=archive/row['blob'];require(blob.is_file() and sha(blob)==row['compressed_sha256'],'missing/changed retained blob:'+name)
        require(decoded(blob)==(row['sha256'],row['bytes']),'lossy retained payload:'+name)
        if file.exists():require(sha(file)==row['sha256'],'raw payload changed before prune:'+name)


def prune(root,archive,index,proof):
    # Validate the entire group before the first unlink: a partial/foreign
    # archive cannot silently erase a later member of the group.
    validate(root,archive,index,proof)
    for name in index['files']:
        file=root/name
        if file.exists():file.unlink()


def restore(root,archive,index):
    for name,row in index['files'].items():
        file=root/name
        require(file.is_relative_to(root) and '..' not in Path(name).parts,'foreign restore path')
        if file.exists():require(sha(file)==row['sha256'],'refuse to replace changed payload');continue
        require(sha(archive/row['blob'])==row['compressed_sha256'],'changed compressed payload')
        temporary=file.with_name(file.name+'.ct12-restore')
        require(not temporary.exists(),'stale restore staging file')
        require(decoded(archive/row['blob'],temporary)==(row['sha256'],row['bytes']),'restore content mismatch')
        os.replace(temporary,file)


def retain(root,archive,eligible,proof):
    require(not archive.exists(),'immutable retention destination already exists')
    archive.mkdir();(archive/'blobs').mkdir();start=time.process_time();wall=time.monotonic();files={};blobs={}
    raw_before=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    for file in sorted(set(eligible)):
        file=file.resolve();require(file.is_relative_to(root.resolve()),'historical/foreign payload cannot be pruned')
        digest=sha(file);blob=archive/'blobs'/(digest+'.gz')
        if not blob.exists():
            with file.open('rb') as src,blob.open('xb') as out:
                with gzip.GzipFile(filename='',fileobj=out,mode='wb',compresslevel=3,mtime=0) as z:shutil.copyfileobj(src,z,1024*1024)
            blobs[digest]=dict(bytes=blob.stat().st_size,sha256=sha(blob))
        files[str(file.relative_to(root))]=dict(sha256=digest,bytes=file.stat().st_size,
            blob=str(blob.relative_to(archive)),compressed_sha256=sha(blob))
    index=dict(schema='LosslessSegmentCASV1',binding=proof['binding'],files=files)
    validate(root,archive,index,proof);neg={}
    for name,mutate in [('missing_blob',lambda q:next(iter(q['files'].values())).update(blob='blobs/absent.gz')),
                        ('changed_plaintext_hash',lambda q:next(iter(q['files'].values())).update(sha256='foreign')),
                        ('foreign_native_binding',lambda q:q.update(binding='foreign'))]:
        q=copy.deepcopy(index);mutate(q)
        try:prune(root,archive,q,proof)
        except (ValueError,RuntimeError,FileNotFoundError) as error:neg[name]=dict(status='pass',reason=str(error))
        else:raise ValueError('retention negative accepted:'+name)
        require(all((root/k).is_file() for k in index['files']),'partial prune after failed validation')
    write_json(archive/'index.json',index);write_json(archive/'required-evidence.json',proof)
    prune(root,archive,index,proof)
    retained=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())+sum(p.stat().st_size for p in archive.rglob('*') if p.is_file())
    comp_cpu=time.process_time()-start;comp_wall=time.monotonic()-wall
    tick=time.process_time();restore(root,archive,index);restore_cpu=time.process_time()-tick
    require(all(sha(root/k)==v['sha256'] for k,v in files.items()),'actual restored bytes differ')
    # Keep the new demonstration's original fields accessible for future audit;
    # the qualified low-retention state was measured before this verified restore.
    return dict(schema='CT12RetentionQualificationV1',status='pass',policy='LosslessSegmentCASV1',
        index_sha256=sha(archive/'index.json'),required_evidence_sha256=sha(archive/'required-evidence.json'),
        original_payload_files=len(files),unique_blobs=len(blobs),raw_segment_bytes=raw_before,
        eligible_raw_bytes=sum(v['bytes'] for v in files.values()),compressed_unique_bytes=sum(v['bytes'] for v in blobs.values()),
        measured_retained_bytes=retained,peak_raw_plus_archive_bytes=raw_before+sum(p.stat().st_size for p in archive.rglob('*') if p.is_file()),
        compression_validation_prune_CPU=comp_cpu,compression_validation_prune_wall=comp_wall,restore_CPU=restore_cpu,
        negative_controls=neg,all_restored_bytes_equal=True,
        final_storage_scope='demonstration restored new payloads after measuring bounded-retention state; no historic archive touched')
