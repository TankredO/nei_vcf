''' setup
'''

import re
import io
from distutils.command.build_ext import build_ext as build_ext_orig
from setuptools import setup, find_packages, Extension

# source: https://stackoverflow.com/a/39671214
__version__ = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open('nei_vcf/__init__.py', encoding='utf_8_sig').read()
).group(1)

# ==== ctypes extensions
class CTypesExtension(Extension):
    '''CTypesExtension'''

class build_ext(build_ext_orig):
    '''build_ext'''
    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypesExtension)
        return super().build_extension(ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return super().get_export_symbols(ext)

    def get_ext_filename(self, ext_name):
        if self._ctypes:
            return ext_name + '.so'
        return super().get_ext_filename(ext_name)

nei_module = CTypesExtension(
    'nei_vcf.lib.nei',
    sources=['nei_vcf/src/nei.cpp'],
    language='c++',
)

vcf_module = CTypesExtension(
    'nei_vcf.lib.vcf',
    sources=['nei_vcf/src/vcf.cpp'],
    language='c++',
)

ext_modules = [
    nei_module,
    vcf_module,
]

install_requires = [
    'numpy',
]

# ====
description = 'Nei (SNP) distance calculation for VCF data.'

long_description = io.open('README.md').read()
long_description_content_type = 'text/markdown'
# ====
setup(
    name='nei_vcf',
    version=__version__,
    packages=find_packages(),
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    author='Tankred Ott',
    platforms=['any'],
    python_requires='>=3.6',
    install_requires=install_requires,
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    # url='',
    entry_points = {
        'console_scripts': [
            'nei_vcf=nei_vcf.commandline:main'
        ],
    },
)
