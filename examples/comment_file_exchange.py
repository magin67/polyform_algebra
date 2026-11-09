"""Synthetic 1C-style comments: JSON snapshot -> activity counts -> result JSON.

This demonstrates transport and stable IDs, not graph centrality or a live 1C connection.
"""
from collections import Counter
from pathlib import Path
import tempfile

from polyform_algebra import get_version
from polyform_algebra.applications.comments import snapshot_digest
from polyform_algebra.connectors.files import (
    CommentSnapshotFileSource, CommentResultFileSink, read_comment_result,
)


def activity_result(snapshot):
    """Count each unique comment once; retain employees with no comments."""
    counts = Counter(row['employee_id'] for row in snapshot['comments'])
    return {
        'schema_version':'1','kind':'comment_results','run_id':'activity-demo-1',
        'source_system':snapshot['source_system'],'snapshot_id':snapshot['snapshot_id'],
        'snapshot_digest':snapshot_digest(snapshot),
        'model':{'id':'comment_count','version':'1','parameters':{}},
        'library_version':get_version(),'created_at':snapshot['exported_at'],
        'status':'completed','metrics':[
            {'employee_id':e['id'],'metric':'comment_count','status':'exact',
             'value':{'type':'rational','numerator':str(counts[e['id']]),'denominator':'1'}}
            for e in snapshot['employees']], 'edges':[], 'messages':[],
    }


def main():
    source = CommentSnapshotFileSource(Path(__file__).parent/'data/comments_snapshot.json')
    snapshot = source.read_snapshot({})
    with tempfile.TemporaryDirectory(prefix='pmg-comments-') as directory:
        sink = CommentResultFileSink(directory,snapshot)
        result = activity_result(snapshot)
        receipt = sink.publish(result,idempotency_key=result['run_id'])
        restored = read_comment_result(receipt['path'],snapshot)
        assert [r['value']['numerator'] for r in restored['metrics']] == ['2','1','0']
        assert sink.publish(result,idempotency_key=result['run_id'])['status'] == 'existing'
        print('File exchange verified: employees 001/002/003 have 2/1/0 comments; retry creates no duplicate.')


if __name__ == '__main__':
    main()
