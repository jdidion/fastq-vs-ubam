from concurrent import futures
import itertools
import gzip
import os
import threading
from typing import Tuple

from ngsindex.utils import BinReader

def nc_record(record):
    A = C = G = T = 0
    for nuc in record[1]:
        if nuc == "A":
            A += 1
        elif nuc == "C":
            C += 1
        elif nuc == "G":
            G += 1
        elif nuc == "T":
            T += 1
        else:
            raise Exception(f"unsupported nucleotide {nuc}")
    
    return (A, C, G, T)


def nc_records(records):
    A = C = G = T = 0
    for record in records:
        for nuc in record[1]:
            if nuc == "A":
                A += 1
            elif nuc == "C":
                C += 1
            elif nuc == "G":
                G += 1
            elif nuc == "T":
                T += 1
            else:
                raise Exception(f"unsupported nucleotide {nuc}")
    
    return (A, C, G, T)


def nc_bytes(b, name_length, seq_length, batch_length, batch_size) -> Tuple[int, int, int, int]:
    A = C = G = T = 0
    for start in range(name_length + 1, batch_length, batch_size):
        for nuc in b[start:(start + seq_length)]:
            if nuc == b'A':
                A += 1
            elif nuc == b'C':
                C += 1
            elif nuc == b'G':
                G += 1
            elif nuc == b'T':
                T += 1
            else:  # the sequence might be whitespace-padded
                break
    
    return (A, C, G, T)


class ThreadSafeNucCounter:
    _lock = threading.Lock()
    _counts = [0, 0, 0, 0]
    
    @staticmethod
    def update(A, C, G, T):
        with ThreadSafeNucCounter._lock:
            ThreadSafeNucCounter._counts[0] += A
            ThreadSafeNucCounter._counts[1] += C
            ThreadSafeNucCounter._counts[2] += G
            ThreadSafeNucCounter._counts[3] += T
    
    @staticmethod
    def get_counts() -> Tuple[int]:
        with ThreadSafeNucCounter._lock:
            return tuple(ThreadSafeNucCounter._counts)


class NcBytesThread(threading.Thread):
    def __init__(
        self, gzfile, read_lock, counter, name_length, seq_length, batch_length, batch_size
    ):
        self.gzfile = gzfile
        self.read_lock = read_lock
        self.counter = counter
        self.name_length = name_length
        self.seq_length = seq_length
        self.batch_length = batch_length
        self.batch_size = batch_size
    
    def run(self):
        while True:
            with self.read_lock:
                b = self.gzfile.read(self.batch_length)
        
            if b:
                self.counter.update(*nc_bytes(
                    b, self.name_length, self.seq_length, self.batch_length, self.batch_size
                ))
            else:
                return


class NcRecordThread(threading.Thread):
    def __init__(self, gzfile, read_lock, counter):
        self.gzfile = gzfile
        self.read_lock = read_lock
        self.counter = counter
    
    def run(self):
        record_itr = zip(*([self.gzfile] * 4))
        while True:
            with self.read_lock:
                try:
                    record = next(record_itr)
                except StopIteration:
                    return

            self.counter.update(*nc_record(record))


class NcRecordsThread(threading.Thread):
    def __init__(self, gzfile, read_lock, counter, batch_size):
        self.gzfile = gzfile
        self.read_lock = read_lock
        self.counter = counter
        self.batch_size = batch_size
    
    def run(self):
        record_itr = zip(*([self.gzfile] * 4))
        while True:
            with self.read_lock:
                chunk = tuple(itertools.islice(record_itr, self.batch_size))
            
            if chunk:
                self.counter.update(*nc_records(chunk))
            else:
                return
        

def nc_fastq_threads(
    fqgz_file, record_lengths=None, batch_size=1, max_workers=None
) -> Tuple[int, int, int, int]:
    gzfile = gzip.open(fqgz_file)
    read_lock = threading.Lock()
    counter = ThreadSafeNucCounter()

    try:
        if record_lengths:
            name_length, seq_length = record_lengths
            record_length = name_length + (2 * seq_length) + 5
            batch_length = record_length * batch_size
            threads = [
                NcBytesThread(
                    gzfile, read_lock, counter, name_length, seq_length, batch_length, batch_size
                )
                for _ in range(max_workers or os.cpu_count())
            ]
        elif batch_size == 1:
            threads = [
                NcRecordThread(gzfile, read_lock, counter)
                for _ in range(max_workers or os.cpu_count())
            ]
        else:
            threads = [
                NcRecordsThread(gzfile, read_lock, counter, batch_size)
                for _ in range(max_workers or os.cpu_count())
            ]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
    finally:
        try:
            gzfile.close()
        except:
            print("error closing fqgz file")
    
    return counter.get_counts()


def nc_fastq_processes(
    fqgz_file, record_lengths=None, batch_size=1, max_workers=None
) -> Tuple[int, int, int, int]:
    pool = futures.ProcessPoolExecutor(max_workers)
    with gzip.open(fqgz_file) as inp:
        if record_lengths:
            name_length, seq_length = record_lengths
            record_length = name_length + (2 * seq_length) + 5
            batch_length = record_length * batch_size
            def submit_next():
                b = inp.read(batch_length)
                if b:
                    return pool.submit(
                        nc_bytes, b, name_length, seq_length, record_length, batch_size
                    )
                else:
                    return None
        else:
            record_itr = zip(*([inp] * 4))
            if batch_size == 1:
                def submit_next():
                    try:
                        return pool.submit(nc_record, next(record_itr))
                    except:
                        return None
            else:
                def submit_next():
                    chunk = tuple(itertools.islice(record_itr, batch_size))
                    if chunk:
                        return pool.submit(nc_records, chunk)
                    else:
                        return None

        fs = []
        done = False
        A = C = G = T = 0
        EXTRA_BATCHES = 1

        def handle_result(f):
            result = f.result()
            nonlocal A, C, G, T
            A += result[0]
            C += result[1]
            G += result[2]
            T += result[3]

        # submit as many batches as workers, plus some extra so the pool is never idle
        for _ in range(max_workers + EXTRA_BATCHES):
            f = submit_next()
            if f:
                fs.extend(f)
            else:
                done = True
                break
        
        # wait for a batch to finish before submitting the next, so we don't use more memory 
        # than we need
        while True:
            if done:
                for f in futures.as_completed(fs):
                    handle_result(f)
                break
            else:
                fs_new = futures.wait(fs, return_when=futures.FIRST_COMPLETED)
                handle_result(fs_new.done[0])
                next_f = submit_next()
                if next_f:
                    fs = fs_new.not_done + (next_f,)
                else:
                    done = True
                    fs = fs_new.not_done
    
    return (A, C, G, T)


def nc_bam(bam, offset) -> Tuple[int, int, int, int]:
    pass


def iter_bgzf_offsets(bam_fileobj):
    reader = BinReader(fileobj=bam_fileobj)
    print(reader.read_vector([("char", 4)]))
    header_len = reader.read_int()
    print(header_len)
    offset = header_len + 12
    print(offset)
    reader.skip(header_len)
    num_ref = reader.read_int()
    print(num_ref)
    for _ in range(num_ref):
        name_len = reader.read_int()
        offset += name_len + 8
    print(offset)
    # read just the BSIZE from each block to compute the next offset
    reader.seek(offset)
    while True:
        yield offset
        reader.skip(16)
        try:
            bsize = reader.read_short()
        except:
            # We're beyond the end of the file
            return
        offset += bsize + 1
        reader.seek(offset)


def nc_bam_processes(
    bam_file, index_file=None, max_workers=None
) -> Tuple[int, int, int, int]:
    pool = futures.ProcessPoolExecutor(max_workers)

    if index_file:
        with open(bam_file, "rb") as bam, open(index_file, "r") as idx:
            for l in idx.readlines():
                pool.submit(nc_bam, bam, int(l.strip()))
    else:
        # TODO: limit the read-ahead to improve cache locality
        with open(bam_file, "rb") as bam:
            for offset in iter_bgzf_offsets(bam):
                pool.submit(nc_bam, bam, offset)
