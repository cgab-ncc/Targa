from setuptools import setup

setup(
    name = 'targa',
    packages = ['targa'], # this must be the same as the name above
    version = '0.0.1',
    description = 'A Python package for identifying targetable genes in cancer with GSEA and KEGG pathway enrichment analysis',
    author = 'Andy Jinseok Lee',
    author_email = 'jinseok.lee@ncc.re.kr',
    url = 'https://github.com/cgab-ncc/Targa',   # use the URL to the github repo
    keywords = ['targetable', 'pathway', 'enrichment', 'genes', 'multiple gene lists', 'cancer', 'genomics'], # arbitrary keywords
    classifiers = [],
    long_description=open('README.md').read(),
    include_package_data=True,
    install_requires=["numpy >= 1.4.0", "scipy >= 0.9.0", "pandas >= 0.17.0", "gseapy >= 0.9.3"],
)