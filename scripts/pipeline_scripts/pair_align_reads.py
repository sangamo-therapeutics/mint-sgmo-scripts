"""Functions and __main__ script to pair R1 and R2, align paired reads to genome using Bowtie2 and filtering for
    MAPQ of 23, make bam file, sort, index

    1. Looks for R1 and R2 .fastq.gz files in "trimmed_fastq" directory and with only string difference
        of "_R1_" and "_R2_" in filename
    2. Makes "aligned" output directory with sub-folder "paired_reads"
    3. Seqkit pairs reads, dropping any unpaired that were individually removed from separate read trimming steps
    4. Performs paired alignment with paired R1 and R2 files, writing sam files to output directory
    5. Use samtools to convert sam to bam files and sort, create index
"""

import os
import re
import subprocess
from pathlib import Path
from typing import Union

THREADS = int(os.getenv("CPUS")) or 16
BT_INDEX_PATH = os.getenv("BOWTIE_INDEX_PATH", "/seqdb/bowtie2")
BT_INDEX = os.getenv("BOWTIE_INDEX", "hg38")
_bt_path = Path(BT_INDEX_PATH) / BT_INDEX
BT_PATH = str(_bt_path)


############################################

def pair_and_align(input_directory: Union[str, Path], output_directory: Union[str, Path, None] = None,
                   index_path: Union[str, Path, None] = BT_PATH, **kwargs):
    """Runs pairing and aligning for all files in

    Args:
        input_directory: Directory containing .fastq.gz
        output_directory: Optional directory for aligned reads; default with be parallel to input dir
        index_path: path to directory and bt2 index files --
        **kwargs:

    """
    input_directory = Path(input_directory)
    output_directory = output_directory or input_directory.parent / "aligned"
    _create_new_dir(output_directory)

    print(f"{output_directory} exists={output_directory.exists()}")
    paired_dir = output_directory / "paired_reads"

    print(f"{paired_dir} exists={paired_dir.exists()}")
    _create_new_dir(paired_dir)
    prior_files = paired_dir.glob("*.*")
    prior_files = set(prior_files)
    #  mode=0o777) ??

    for r1, r2 in _find_readpairs(input_directory):
        _make_paired_fastq(r1, r2, paired_dir, **kwargs)
        paired_r1 = paired_dir / r1.name
        paired_r2 = paired_dir / r2.name
        # paired_r1.chmod(0o777)  # necessary since container may assign to root
        # paired_r2.chmod(0o777)
        _run_aligner(paired_r1, paired_r2, index_path, output_directory=output_directory)
    # output_files = output_directory.glob("*.*")
    # for f in output_files:
    #     if f not in prior_files:
    #         for f in new_files:
    #             try:
    #                 os.chown(f, HOST_UID, HOST_GID)
    #             except Exception as e:
    #                 print(f"unable to change permissions on file: {f} {e}")


def _create_new_dir(output_directory):
    if not output_directory.exists():
        print(f"############## Creating output directory:{output_directory} #########")
        output_directory.mkdir(mode=0o777, exist_ok=True)
        if HOST_UID:
            print(f"Setting {output_directory} to user {HOST_UID}")
            try:
                os.chown(output_directory, HOST_UID, HOST_GID)
            except Exception as e:
                print("Unable to set user.")


def convert_sam_to_bam(sam_input: Union[str, Path], bam_output: Union[str, Path, None] = None):
    """Converts sam to bam and sorts and indexes bam

    Args:
        sam_input: full path to sam file
        bam_output: optional full path to an output bam file.
            If not passed, the bam will be written in the same directory with only the extension changed

    """
    sam_input = Path(sam_input)
    bam_file = bam_output or sam_input.parent / sam_input.name.replace(".sam", ".bam")
    bam_file_sorted_name, bam_file_sorted_coord = _get_sorted_bam_paths(sam_input)
    with open(bam_file_sorted_name, "wb") as output_file:
        process1 = subprocess.Popen(
            ["samtools", "view", "-@", "16", "-bh", "-q", "23", str(sam_input)],
            stdout=subprocess.PIPE)
        process2 = subprocess.Popen(["samtools", "sort", "-@", "16", "-n", "-o", "-", "-"], stdin=process1.stdout,
                                    stdout=output_file)
        process1.stdout.close()
        process2.communicate()

    with open(bam_file_sorted_coord, "wb") as output_file:
        process3 = subprocess.Popen(["samtools", "sort", "-@", "16", "-o", "-", bam_file_sorted_name],
                                    stdout=output_file)
        process3.communicate()

        # Index the BAM file sorted by genomic coordinates
    subprocess.run(["samtools", "index", "-@", "16", bam_file_sorted_coord])


def _find_readpairs(input_directory: Union[str, Path]):
    """ Searches a directory and yields pairs of matching reads by readname
    Names must be identical other than _R1_ _R2_ designation

    Files must be in same directory

    Args:
        input_directory: path to location of the fastq.gz files

    Yields (tuple): paired R1 / R2 files

    """
    input_directory = Path(input_directory)
    read2_files = input_directory.glob('*_R2_*.fastq.gz')
    # read2_files = [f for f in infiles if "_R2_" in f]
    for r2 in read2_files:
        r1 = r2.parent / r2.name.replace("_R2_", "_R1_")
        if r1.exists():
            yield r1, r2
        else:
            print(f"{r2.name} :file is missing read1 file")


def _make_paired_fastq(read1: Union[str, Path], read2: Union[str, Path], output_directory: Union[str, Path],
                       threads: int = THREADS, **kwargs):
    """For a pair of read files, pairs them with seqkit and writes new sorted paired reads

    Args:
        read1: path to R1 file
        read2: path to R2 file
        output_directory: path or Path to location to write paired fastq files
        threads: number of cpus to utilize

    """
    output_dir = Path(output_directory)
    output_dir.mkdir(exist_ok=True)
    seqkit_command = ["seqkit", "pair",
                      "-1", str(read1),
                      "-2", str(read2),
                      "-O", str(output_dir),
                      "-j", f"{threads}"]
    try:
        subprocess.run(seqkit_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit pair: {e}")
        return


def _get_reads_base(read_path: Union[str, Path], ):
    """Gets the 'basename of the fastq.gz file, removing the read number """
    read_path = Path(read_path)
    base = re.sub("_R\d_", "_", read_path.name)
    if ".fastq.gz" in base:
        base = base.replace(".fastq.gz", "")
    else:
        base = base.split(".")[0]
    return base


def _run_aligner(read1_path: Union[str, Path], read2_path: Union[str, Path], index: str,
                 output_directory: Union[str, Path]):
    """Executes bowtie2 aligner on two files and writes resulting ordered reads to an output directory """
    base_name = _get_reads_base(read1_path)
    output_directory = Path(output_directory)
    output_directory.mkdir(exist_ok=True)
    output_sam = output_directory / f"{base_name}.sam"
    old_files = output_directory.glob(f"{base_name}*")
    old_files = list(old_files)
    try:
        # print(f"checking index at {index}")
        assert Path(index).parent.exists()
    except AssertionError:
        print(f"I could not find path to {index}")
        raise
    bowtie2_command = ["bowtie2", "-x", index, "-1", str(read1_path), "-2", str(read2_path),
                       "--very-sensitive-local", "--threads", "16", "--no-unal", "-S",
                       str(output_sam)]
    print(" ".join(bowtie2_command))
    qc_file = output_directory / f"{base_name}_bowtieQC.txt"
    with open(qc_file, "a") as f:
        subprocess.run(bowtie2_command, stdout=subprocess.PIPE, stderr=f)
    old_files = output_directory.glob(f"{base_name}*")
    old_files = list(old_files)


def _get_sorted_bam_paths(in_path: Union[str, Path]):
    """Create standard paths and filenames from an input sam file

    Args:
        in_path: path to source sam or bam file

    Returns (tuple): Path objects for the _sorted_name and _sorted_coord bam files

    """
    in_path = Path(in_path)
    base_name = in_path.name.replace(in_path.suffix, "")
    bam_file_sorted_name = in_path.parent / f"{base_name}_sorted_name.bam"
    bam_file_sorted_coord = in_path.parent / f"{base_name}_sorted_coord.bam"
    return bam_file_sorted_name, bam_file_sorted_coord


def _set_newfile_permissions(directory, prior_files=None, host_uid=None, host_gid=None):
    if not host_uid:
        return
    prior_files = prior_files or []
    new_files = directory.glob(f"*.*")
    new_files = [f for f in new_files if f not in prior_files]
    for nf in new_files:
        try:
            os.chown(nf, host_uid, host_gid)
        except Exception as e:
            print(f"unable to change permissions on file: {nf} {e}")


if __name__ == "__main__":
    try:
        HOST_UID = int(os.getenv("HOST_UID"))
    except (TypeError, ValueError):
        HOST_UID = None
    try:
        HOST_GID = int(os.getenv("HOST_GID"))
    except (TypeError, ValueError):
        HOST_GID = None
    input_dirname = os.getenv("DEFAULT_PIPELINE_IN", "/fastq_reads")
    output_dirname = os.getenv("DEFAULT_PIPELINE_OUT", "/aligned")
    input_dir = Path(input_dirname)
    output_dir = Path(output_dirname)

    try:
        assert input_dir.exists()
    except AssertionError:
        print(f"Directory for input {str(input_dir)} does not exist.")
        raise
    print("Executing...")
    try:
        assert output_dir.exists()
    except AssertionError:
        print(f"Directory for output {str(output_dir)} does not exist.")
        raise
    print("Executing...")

    index_path = BT_PATH
    try:
        print(f"INDEX_PATH = {index_path}")
        assert Path(index_path).parent.exists()
    except AssertionError:
        print(f"something is wrong with the index path {index_path}")
        raise

    # ############ run pairing and aligning
    paired_dir = output_dir / "paired_reads"
    prior_pairfiles = paired_dir.glob("*.*")
    prior_pairfiles = set(prior_pairfiles)

    pair_and_align(input_dir, output_directory=output_dir, index_path=index_path,
                   index_name=BT_INDEX)  # using default output dir and index path and threads

    _set_newfile_permissions(paired_dir, prior_files=prior_pairfiles, host_uid=HOST_UID, host_gid=HOST_GID)

    # ######## sam files should now be in the "output_dir"
    samfiles = output_dir.glob("*.sam")
    samfiles = list(samfiles)
    if not samfiles:
        print("No sam files were found for conversion.")
    prior_files = output_dir.glob("*.*")
    prior_files = set(prior_files)
    for sam_file in samfiles:
        # sam_file.chmod(0o666)
        basename = _get_reads_base(sam_file)
        try:
            os.chown(sam_file, HOST_UID, HOST_GID)
        except Exception as e:
            print(f"unable to change permissions on file: {sam_file} {e}")
        convert_sam_to_bam(sam_file)
        # new_files = output_dir.glob(f"{basename}*")
        # new_files = [f for f in new_files if f not in prior_files]
        # for f in new_files:
        #     try:
        #         os.chown(f, HOST_UID, HOST_GID)
        #     except Exception as e:
        #         print(f"unable to change permissions on file: {f} {e}")
    _set_newfile_permissions(output_dir, prior_files, host_uid=HOST_UID, host_gid=HOST_GID)
