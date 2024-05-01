from setuptools import setup

setup(
  name = 'NODE_deconvolution',         
  packages = ['NODE_deconvolution'],  
  version = '0.1.3',      
  license='MIT',        
  description = 'A Python program for deconvolution spatial transcriptomics data and inference of spatial communication',   
  author = 'Zedong Wang',                   
  author_email = 'wangzedong23@mails.ucas.ac.cn',     
  url = 'https://github.com/wzdrgi/NODE',   
  keywords = ["spatial transcriptomics", 'deconvolution ', 'spatial communication'],   
  install_requires=[           
          'scipy',  
          'numpy',
          'tqdm',
          'scikit-opt'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',   
    'Programming Language :: Python :: 3.5', 
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3.12',
  ],
)
