from distutils.core import setup, Extension 
from Cython.Distutils import build_ext
from Cython.Build import cythonize 

ext_modules = [
               Extension("traverse_contours", 
                         ["traverse_contours.pyx"],
                         language="c",
                         include_dirs = ["/net/eichler/vol6/software/modules-sw/anaconda/2.1.0/Linux/RHEL6/x86_64/envs/python2/lib/python2.7/site-packages/numpy/core/include/"]
                         ),
               Extension("get_windowed_variance", 
                         ["get_windowed_variance.pyx"],
                         language="c",
                         include_dirs = ["/net/eichler/vol6/software/modules-sw/anaconda/2.1.0/Linux/RHEL6/x86_64/envs/python2/lib/python2.7/site-packages/numpy/core/include/"]
                         ),
               Extension("c_hierarchical_edge_merge",
                         ['c_hierarchical_edge_merge.pyx',
                         './heap/heap.cc',
                         './edge.cc',
                         './float_heap.cc'],
                        language = "c++",
                        include_dirs = ["/net/eichler/vol6/software/modules-sw/anaconda/2.1.0/Linux/RHEL6/x86_64/envs/python2/include/python2.7/", "heap/"]
                        )
              ]

setup( 
    cmdclass = {'build_ext':build_ext},
    name = 'c_SSF',
    ext_modules = ext_modules 
    )
