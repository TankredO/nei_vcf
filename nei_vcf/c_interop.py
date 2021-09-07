''' C library bindings
'''

import pathlib
import ctypes
from typing import Tuple
import numpy as np
import numpy.ctypeslib as ctl


# load the shared libraries
vcf_libfile = pathlib.Path(__file__).parent / 'lib' / 'vcf.so'
vcf_lib = ctypes.CDLL(str(vcf_libfile))

nei_libfile = pathlib.Path(__file__).parent / 'lib' / 'nei.so'
nei_lib = ctypes.CDLL(str(nei_libfile))


# get VCF file "height" (number of loci) and "width" (number of samples)
get_vcf_hw = vcf_lib.vcf_stats
get_vcf_hw.argtypes = [
    ctypes.c_char_p, # file name
    ctl.ndpointer(np.int32, flags='aligned, c_contiguous'), # shape array (out)
    ctypes.c_int, # n unused columns
]

def get_vcf_shape(
    vcf_file: str,
    n_col_skip: int = 9,
) -> np.ndarray:
    vcf_file = vcf_file.encode('utf-8')
    shape = np.array([0, 0], dtype=np.int32)
    
    get_vcf_hw(vcf_file, shape, n_col_skip)

    return shape

# get variant information from vcf file
get_vcf_variants = vcf_lib.get_variants
get_vcf_variants.argtypes = [
    ctypes.c_char_p, # file name (in)
    ctl.ndpointer(np.double, flags='aligned, c_contiguous'), # variants array (out)
    ctypes.c_char_p, # sample_name_buffer (out)
    ctl.ndpointer(np.int32, flags='aligned, c_contiguous'), # shape array (in)
    ctypes.c_int, # n unused columns
    ctypes.c_int, # entry index
]
MAX_SAMPLE_NAME_LEN = 256

def get_variants(
    vcf_file: str,
    n_col_skip: int = 9,
    entry_index: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    shape = get_vcf_shape(vcf_file, n_col_skip)
    vcf_file = vcf_file.encode('utf-8')

    sample_names_buffer = ctypes.create_string_buffer(b"", size=(MAX_SAMPLE_NAME_LEN + 1) * shape[1])
    variants = np.zeros(shape=shape[0]*shape[1]*4, dtype=np.double)
    get_vcf_variants(
        vcf_file,
        variants,
        sample_names_buffer,
        shape,
        n_col_skip,
        entry_index
    )

    sample_names = np.array([
        s.rstrip() for s in sample_names_buffer.raw.decode('utf-8').split('\x00') if s
    ])
    variant_arr = variants.reshape((shape[0], shape[1], 4))

    return variant_arr, sample_names

# calculate Nei distance based on variants
get_nei_dist = nei_lib.nei_dist
get_nei_dist.argtypes = [
    ctl.ndpointer(np.double, flags='aligned, c_contiguous'), # variants array (in)
    ctl.ndpointer(np.int32, flags='aligned, c_contiguous'), # shape array (in)
    ctl.ndpointer(np.double, flags='aligned, c_contiguous'), # dist array (out)
    ctl.ndpointer(np.int32, flags='aligned, c_contiguous'), # count array (out)
]

def nei(
    variants: np.ndarray,
) -> np.ndarray:
    shape = np.array(variants.shape, dtype=np.int32)
    distances = np.zeros(shape=shape[1] * shape[1], dtype=np.double)
    counts = np.zeros(shape=shape[1] * shape[1], dtype=np.int32)

    get_nei_dist(variants, shape, distances, counts)

    distances = distances.reshape((shape[1], shape[1]))
    counts = counts.reshape((shape[1], shape[1]))

    distances = distances / counts

    return distances

def nei_vcf(
    vcf_file: str,
    n_col_skip: int = 9,
    entry_index: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    variants, sample_names = get_variants(vcf_file, n_col_skip, entry_index)
    distances = nei(variants)

    return distances, sample_names
