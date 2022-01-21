# FASTQ vs UBAM
## Background
## Usability
* Sample- and read-level metadata
* Single file for both single- and paired-end reads
* Different query names for reads 1 and 2 - store suffix(es) in tags
## File size
* FASTQ will always be superior unless BAM allows different compression schemes
  * More compression options for FASTQ
## Performance
### Naive/serial decompression speed
### Parallel decompression/parsing/processing
* Simple task (compute nucleotide composition)
* Add random sleep to simulate more complex task
#### Multi-threading (FASTQ and uBAM)
* Share single file handle across threads
* Synchronize on file to read batch
* Record padding for deferred parsing
#### Multi-processing (FASTQ and uBAM)
* Main process reads file and adds batches to queue
* Sub-processes pull batches off queue and process
* Record padding for deferred parsing (main process vs subprocess)
#### Unsynchronized parallel processing of bgzf blocks (uBAM only)
* Main process adds block offsets and lengths to queue
  * Threads/subprocesses pull block coordinates off the queue, read, and process
* Use of index vs main process reading block offsets at runtime
* Memory mapping
  * Limit read-ahead to keep memory to fixed size
## Discussion
### Alternatives
* BGZIP FASTQ
* New custom format
  * Column-oriented
    * Different compression scheme for each column
    * Decompress columns in parallel
    * Don't decompress unused columns
    * Single file
      * Eliminate redundant storage of query name
      * Separate columns for reads 1 and 2 (self-describing)
    * Optimized storage of long reads
    * Easy interop w/ Arrow