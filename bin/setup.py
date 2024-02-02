from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='tetrimmer',
    version='1.1.5',
    cmdclass=versioneer.get_cmdclass(),
    description="a tool to replace transposable element manual curation. TE Trimmer won't do TE de novo annotation but use the output from other annotation tools like RepeatModeler, REPET, and EDTA",
    license="MIT",
    author="Jiangzhao Qian",
    author_email='jqian@bio1.rwth-aachen.de',
    url='https://github.com/qjiangzhao/TETrimmer',
    packages=['tetrimmer'],
    entry_points={
        'console_scripts': [
            'tetrimmer=tetrimmer.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='TETrimmer',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
