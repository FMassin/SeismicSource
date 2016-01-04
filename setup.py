try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project aim at probabilistic inversion of the source parameters of earthquakes',
    'author': 'Frederick Massin',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'fred.massin@gmail.com',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['NnK'],
    'scripts': [],
    'name': 'Naino-Kami'
}

setup(**config)
