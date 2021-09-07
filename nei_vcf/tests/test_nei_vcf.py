import pathlib
import shutil
import tempfile

import numpy as np
import pytest

from nei_vcf import get_variants, nei_vcf, nei
from nei_vcf import commandline

def test_simple():
    vcf_file = str(pathlib.Path(__file__).parent / 'test_data.vcf')

    variants, sample_names = get_variants(vcf_file)
    print('variants:')
    print(variants)
    print('sample_names:')
    print(sample_names)

    distances = nei(variants)
    print('distances:')
    print(distances)

    distances2, sample_names2 = nei_vcf(vcf_file)
    print('distances2:')
    print(distances2)
    print('sample_names2:')
    print(sample_names2)

    assert np.abs(distances - distances2).sum() < 1e-8
    assert (np.array(sample_names) == np.array(sample_names2)).all()



def test_commandline(capsys):
    with pytest.raises(SystemExit) as wrapped_err:
        commandline.main(['-h'])
    out, err = capsys.readouterr()
    assert wrapped_err.value.code == 0

    with pytest.raises(SystemExit) as wrapped_err:
        commandline.main(['-v'])
    assert wrapped_err.value.code == 0

    vcf_file = str(pathlib.Path(__file__).parent / 'test_data.vcf')
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_file = str(pathlib.Path(tmp_dir).joinpath('test_data.dist'))
        commandline.main([vcf_file, '-o', out_file])

    with tempfile.TemporaryDirectory() as tmp_dir:
        vcf_file_copy = str(pathlib.Path(tmp_dir).joinpath('test_data.vcf'))
        shutil.copy(vcf_file, vcf_file_copy)
        out_file = str(pathlib.Path(tmp_dir).joinpath('test_data.dist'))
        commandline.main([vcf_file_copy, '-o', out_file])

    