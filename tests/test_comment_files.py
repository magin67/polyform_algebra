"""Applied exchange contract: identity, completeness and atomic result publication."""
import copy
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import pytest
from polyform_algebra.applications.comments import DataContractError, validate_snapshot, validate_result, snapshot_digest
from polyform_algebra.connectors.files import CommentSnapshotFileSource, CommentResultFileSink, PublicationConflict, read_comment_result
from examples.comment_file_exchange import activity_result


@pytest.fixture
def snapshot():
    return CommentSnapshotFileSource(Path(__file__).resolve().parents[1]/'examples/data/comments_snapshot.json').read_snapshot({})


def test_roundtrip_preserves_ids_unicode_and_zero_activity(snapshot,tmp_path):
    """Synthetic input has counts 2,1,0; employee/document ID 001 are distinct typed entities."""
    path=tmp_path/'source.json'
    path.write_text('\ufeff'+json.dumps(snapshot,ensure_ascii=False),encoding='utf-8')
    source=CommentSnapshotFileSource(path)
    loaded=source.read_snapshot({'snapshot_id':snapshot['snapshot_id']})
    assert loaded==snapshot and loaded['employees'][0]['id']=='001'
    loaded['employees'][0]['label']='changed'
    assert source.read_snapshot({})==snapshot
    result=activity_result(snapshot)
    sink=CommentResultFileSink(tmp_path/'results',snapshot)
    receipt=sink.publish(result,idempotency_key=result['run_id'])
    assert read_comment_result(receipt['path'],snapshot)==result
    assert [x['value']['numerator'] for x in result['metrics']]==['2','1','0']
    assert sink.publish(copy.deepcopy(result),idempotency_key=result['run_id'])['status']=='existing'
    assert len(list((tmp_path/'results').iterdir()))==1


@pytest.mark.parametrize('defect', ['incomplete','schema','duplicate_event','duplicate_employee','unknown_employee','unknown_document','numeric_id','naive_time','end_time','bad_period','early_export','unknown_field','null_comments'])
def test_snapshot_rejects_ambiguous_or_incomplete_data(snapshot,defect):
    """Reject malformed input instead of silently dropping records or double counting comments."""
    if defect=='incomplete': snapshot['complete']=False
    elif defect=='schema': snapshot['schema_version']='2'
    elif defect=='duplicate_event': snapshot['comments'].append(snapshot['comments'][0])
    elif defect=='duplicate_employee': snapshot['employees'].append(snapshot['employees'][0])
    elif defect=='unknown_employee': snapshot['comments'][0]['employee_id']='missing'
    elif defect=='unknown_document': snapshot['comments'][0]['document_id']='missing'
    elif defect=='numeric_id': snapshot['employees'][0]['id']=1
    elif defect=='naive_time': snapshot['comments'][0]['occurred_at']='2026-10-15T12:00:00'
    elif defect=='end_time': snapshot['comments'][0]['occurred_at']=snapshot['period']['end']
    elif defect=='bad_period': snapshot['period']['end']=snapshot['period']['start']
    elif defect=='early_export': snapshot['exported_at']=snapshot['period']['start']
    elif defect=='unknown_field': snapshot['comments'][0]['deleted']=True
    else: snapshot['comments']=None
    with pytest.raises(DataContractError): validate_snapshot(snapshot)


def test_empty_snapshot_and_offsets(snapshot):
    """Complete empty periods are valid; timestamps compare as instants across UTC offsets."""
    snapshot['comments'][0]['occurred_at']='2026-09-30T20:00:00Z'
    assert validate_snapshot(snapshot)
    snapshot['comments']=[]
    assert validate_snapshot(snapshot)['employees']==snapshot['employees']
    snapshot['employees']=[]; snapshot['documents']=[]
    assert validate_snapshot(snapshot)['comments']==[]


@pytest.mark.parametrize('raw', [b'{',b'{"a":1,"a":2}',b'{"x":NaN}',b'{"x":Infinity}',b'\xff'])
def test_bad_json_is_not_an_empty_snapshot(tmp_path,raw):
    """Broken encodings, duplicate keys and nonstandard numbers fail before processing."""
    path=tmp_path/'bad.json'; path.write_bytes(raw)
    with pytest.raises(DataContractError): CommentSnapshotFileSource(path).read_snapshot({})


def test_file_limits_query_and_missing_file(snapshot,tmp_path):
    """No implicit filter, truncation or successful empty response on a missing file."""
    path=tmp_path/'s.json'; path.write_text(json.dumps(snapshot),encoding='utf-8')
    with pytest.raises(DataContractError): CommentSnapshotFileSource(path,max_bytes=10).read_snapshot({})
    for query in ({'snapshot_id':'wrong'},{'department':'x'}):
        with pytest.raises(DataContractError): CommentSnapshotFileSource(path).read_snapshot(query)
    with pytest.raises(FileNotFoundError): CommentSnapshotFileSource(tmp_path/'missing').read_snapshot({})


@pytest.mark.parametrize('defect',['digest','source','snapshot','employee','duplicate','status','float_exact','zero_denominator','undefined','negative_edge','loop','reverse_duplicate'])
def test_result_binding_and_value_validation(snapshot,defect):
    """A score must retain source identity, explicit precision and valid employee references."""
    result=activity_result(snapshot)
    if defect=='digest': result['snapshot_digest']='bad'
    elif defect=='source': result['source_system']='other'
    elif defect=='snapshot': result['snapshot_id']='other'
    elif defect=='employee': result['metrics'][0]['employee_id']='missing'
    elif defect=='duplicate': result['metrics'].append(result['metrics'][0])
    elif defect=='status': result['status']='partial'
    elif defect=='float_exact': result['metrics'][0]['value']={'type':'number','value':0.5}
    elif defect=='zero_denominator': result['metrics'][0]['value']['denominator']='0'
    elif defect=='undefined': result['metrics'][0].update(status='undefined',value=None)
    else:
        result['edges']=[{'source_id':'001','target_id':'001' if defect=='loop' else '002','weight':{'type':'rational','numerator':'-1' if defect=='negative_edge' else '1','denominator':'2'}}]
        if defect=='reverse_duplicate': result['edges'].append({'source_id':'002','target_id':'001','weight':{'type':'number','value':0.5}})
    with pytest.raises(DataContractError): validate_result(result,snapshot)


def test_undefined_and_infinity_are_not_zero(snapshot):
    """Undefined metrics require a reason; infinity is typed rather than nonstandard JSON."""
    r=activity_result(snapshot)
    r['metrics'][0].update(status='undefined',value=None,reason='isolated')
    r['metrics'][1]['value']={'type':'positive_infinity'}
    r['metrics'][2].update(status='approximate',value={'type':'number','value':0.25})
    assert validate_result(r,snapshot)==r


def test_conflict_and_path_safety(snapshot,tmp_path):
    """A conflicting retry cannot overwrite the first result; external IDs cannot form paths."""
    result=activity_result(snapshot); result['run_id']='../../outside'
    sink=CommentResultFileSink(tmp_path,snapshot)
    receipt=sink.publish(result,idempotency_key=result['run_id'])
    target=Path(receipt['path']); original=target.read_bytes()
    assert target.parent==tmp_path
    result['metrics'][0]['value']['numerator']='99'
    with pytest.raises(PublicationConflict): sink.publish(result,idempotency_key=result['run_id'])
    assert target.read_bytes()==original
    with pytest.raises(DataContractError): sink.publish(result,idempotency_key='different')
    target.write_text('{broken',encoding='utf-8')
    with pytest.raises(PublicationConflict): sink.publish(result,idempotency_key=result['run_id'])


def test_concurrent_publication(snapshot,tmp_path):
    """Concurrent equal retries produce one complete file, with exactly one creator."""
    result=activity_result(snapshot); sink=CommentResultFileSink(tmp_path,snapshot)
    with ThreadPoolExecutor(max_workers=4) as pool:
        receipts=list(pool.map(lambda _:sink.publish(result,idempotency_key=result['run_id']),range(8)))
    assert sum(r['status']=='created' for r in receipts)==1
    assert len(list(tmp_path.iterdir()))==1
    assert read_comment_result(receipts[0]['path'],snapshot)==result


def test_failed_publication_cleans_temporary_files(snapshot,tmp_path,monkeypatch):
    """A filesystem publication failure exposes no final result and leaves no partial temp file."""
    import polyform_algebra.connectors.files as files
    def fail(*args): raise OSError('simulated publication failure')
    monkeypatch.setattr(files.os,'link',fail)
    result=activity_result(snapshot)
    with pytest.raises(OSError): CommentResultFileSink(tmp_path,snapshot).publish(result,idempotency_key=result['run_id'])
    assert list(tmp_path.iterdir())==[]


def test_snapshot_digest_and_detached_sink(snapshot,tmp_path):
    """Results bind to snapshot content, not just a reused human-readable snapshot ID."""
    old=copy.deepcopy(snapshot); sink=CommentResultFileSink(tmp_path,snapshot)
    snapshot['employees'][0]['label']='new label'
    assert snapshot_digest(old)!=snapshot_digest(snapshot)
    result=activity_result(old)
    assert sink.publish(result,idempotency_key=result['run_id'])['status']=='created'
    with pytest.raises(DataContractError): validate_result(result,snapshot)


def test_limits_readonly_mappings_and_invalid_unicode(snapshot,tmp_path):
    """Mapping inputs remain compatible with connector contracts; byte limits apply before publication."""
    from types import MappingProxyType
    assert validate_snapshot(MappingProxyType(snapshot))==snapshot
    result=activity_result(snapshot)
    with pytest.raises(DataContractError):
        CommentResultFileSink(tmp_path,snapshot,max_bytes=10).publish(result,idempotency_key=result['run_id'])
    assert list(tmp_path.iterdir())==[]
    snapshot['employees'][0]['label']='\ud800'
    with pytest.raises(DataContractError): validate_snapshot(snapshot)
