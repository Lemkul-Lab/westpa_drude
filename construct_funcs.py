#!/usr/bin/env python
import numpy

def construct_auxdata1(n_iter, iter_group):
  dim1 = iter_group['auxdata/dih1']
  dim2 = iter_group['auxdata/dih2']
  dataset = numpy.dstack((dim1, dim2))
  return dataset
