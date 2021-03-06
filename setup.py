from setuptools import setup, find_packages

setup(
    name = 'optimModels',
    version = '0.0.16',
    package_dir = {'': 'src'},
    packages = find_packages('src'),
    include_package_data = True,
    install_requires = ['inspyred', 'framed', 'geckopy', 'cobra', 'optlang', 'pandas', 'numpy'],
    author = 'Sara Correia',
    author_email = 'sarag.correia@gmail.com',
    description = 'optimModels - strain optimization',
    license = 'Apache License Version 2.0',
    keywords = 'strain design',
    url = 'https://github.com/saragcorreia/optimModels.git',
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    classifiers = ['Topic :: Utilities', 'Programming Language :: Python :: 3.6', ],
    )
