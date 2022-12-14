#!/usr/bin/env python3
# encoding: utf-8

from distutils.core import setup, Extension


def configuration(parent_package='', top_path=None):
      #import numpy
      from numpy.distutils.misc_util import Configuration
      from numpy.distutils.misc_util import get_info

      config = Configuration('',
                             parent_package,
                             top_path)
      config.add_extension('check_alignment',
                           ['src/check_alignment.cpp', 'src/models.cpp'],
                           extra_objects=['Dependencies/htslib-1.9/libhts.a'],
                           extra_compile_args=['-std=c++14'], language='c++',
                           extra_link_args=['-lz', '-lpthread'],
                           include_dirs=['include', 'Dependencies/htslib-1.9/htslib'])

      return config

if __name__ == "__main__":
      from numpy.distutils.core import setup

      setup(
        name='check_alignment',
        version='0.0.1',
        configuration=configuration)
