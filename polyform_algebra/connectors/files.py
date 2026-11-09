"""UTF-8 JSON file exchange for the employee/document/comment contract.

Publication uses an atomic, no-replace hard link in the destination directory.
Filesystems without hard-link support fail explicitly instead of risking a
partial file or overwriting an earlier result. No database connection is opened.
"""
import hashlib
from collections.abc import Mapping
import json
import os
from pathlib import Path
import tempfile

from .base import DataSource, ResultSink
from ..applications.comments import (
    DataContractError, canonical_json, validate_snapshot, validate_result,
)


class PublicationConflict(DataContractError):
    """A run key already names a different or damaged published result."""


def _read(path, max_bytes):
    with Path(path).open('rb') as stream:
        raw = stream.read(max_bytes+1)
    if len(raw) > max_bytes:
        raise DataContractError('JSON file exceeds configured byte limit')
    def pairs(items):
        result = {}
        for k,v in items:
            if k in result:
                raise DataContractError('duplicate JSON object key')
            result[k] = v
        return result
    def constant(value):
        raise DataContractError('nonfinite JSON number: '+value)
    try:
        return json.loads(raw.decode('utf-8-sig'), object_pairs_hook=pairs, parse_constant=constant)
    except (UnicodeError, ValueError, RecursionError) as error:
        raise DataContractError('invalid UTF-8 JSON file') from error


def _limit(value):
    if type(value) is not int or value <= 0:
        raise DataContractError('max_bytes must be a positive integer')
    return value


class CommentSnapshotFileSource(DataSource):
    """Read one complete snapshot; no implicit filtering or incremental merge."""

    def __init__(self, path, *, max_bytes=16*1024*1024):
        self.path = Path(path)
        self.max_bytes = _limit(max_bytes)

    def read_snapshot(self, query):
        """Accept {} or an expected snapshot_id; return freshly validated detached data."""
        if not isinstance(query, Mapping) or set(query)-{'snapshot_id'}:
            raise DataContractError('query supports only expected snapshot_id')
        snapshot = validate_snapshot(_read(self.path,self.max_bytes))
        if 'snapshot_id' in query and query['snapshot_id'] != snapshot['snapshot_id']:
            raise DataContractError('unexpected snapshot_id')
        return snapshot


class CommentResultFileSink(ResultSink):
    """Publish immutable result files bound to a supplied validated snapshot.

    The idempotency key must equal run_id. Filenames hash (source_system, run_id)
    so retries share one path and external IDs cannot become path traversal.
    Existing files are never replaced. Equal retries return status 'existing'.
    """

    def __init__(self, directory, snapshot, *, max_bytes=16*1024*1024):
        self.directory = Path(directory)
        self._snapshot = validate_snapshot(snapshot)
        self.max_bytes = _limit(max_bytes)

    def publish(self, result, *, idempotency_key):
        value = validate_result(result,self._snapshot)
        if type(idempotency_key) is not str or idempotency_key != value['run_id']:
            raise DataContractError('idempotency_key must equal run_id')
        data = (canonical_json(value)+'\n').encode('utf-8')
        if len(data) > self.max_bytes:
            raise DataContractError('result exceeds configured byte limit')
        key = canonical_json([value['source_system'],value['run_id']]).encode('utf-8')
        target = self.directory / (hashlib.sha256(key).hexdigest()+'.json')
        self.directory.mkdir(parents=True,exist_ok=True)
        status = 'created'
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=self.directory,prefix='.pmg-',suffix='.tmp',delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
            try:
                os.link(temporary,target)
            except FileExistsError:
                try:
                    old = validate_result(_read(target,self.max_bytes),self._snapshot)
                except DataContractError as error:
                    raise PublicationConflict('existing result is invalid or belongs to another snapshot') from error
                if canonical_json(old) != canonical_json(value):
                    raise PublicationConflict('run_id already published with a different result')
                status = 'existing'
        finally:
            if temporary is not None:
                temporary.unlink()
        return {'status':status,'run_id':value['run_id'],'path':str(target.resolve()),
                'sha256':hashlib.sha256(canonical_json(value).encode('utf-8')).hexdigest()}


def read_comment_result(path, snapshot, *, max_bytes=16*1024*1024):
    """Read the same JSON contract returned to 1C and verify its source binding."""
    return validate_result(_read(path,_limit(max_bytes)),snapshot)
