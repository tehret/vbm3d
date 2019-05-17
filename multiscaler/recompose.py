#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>

import os
import sys

import numpy as np
import piio
import pywt


def pick_option(argv, option, default):
   # it's a parameter or just a flag?
   if default == '':
      id = 0
   else:
      id = 1

   for i in range(len(argv)-id):
      if argv[i]  == '-'+option:
         r = argv[i+id]
         argv.pop(i+id)
         if id:
            argv.pop(i)
         return r
   return default


def recompose(prefix, levels, suffix, wtype, cur=0):

   print prefix+str(cur)+suffix
   image = piio.read(prefix+str(cur)+suffix)
   print image.shape

   if levels==1: # if last level
      return image

   LL = 2 * recompose(prefix,levels-1,suffix, wtype, cur+1) ## 2 compensates the factor used in decompose 

   coeffs2  = pywt.dwt2(image,wtype,axes=(0,1))

   if coeffs2[0].shape != LL.shape:
      print "ERROR: Are you sure that \'%s\' is the right wavelet? Sizes don't match!"%wtype
      exit()

   ncoeffs2 = (LL, coeffs2[1])

   ret = pywt.idwt2(ncoeffs2,wtype, axes=(0,1))

   # adjust size of the upsampled image
   ret = ret[:image.shape[0],:image.shape[1],:] 
   print ret.shape
   return ret


def main():
   # verify input
   wtype = pick_option(sys.argv, 't', 'db7')
   if len(sys.argv) > 4:
      prefix = sys.argv[1]
      levels = int(sys.argv[2])
      suffix = sys.argv[3]
      oname  = sys.argv[4]
   else:
      print("Incorrect syntax, use:")
      print("  > " + sys.argv[0] + " prefix levels suffix output [-t type]")
      print("        multiscale recomposition using wavelets")
      print("        available types:")
      for family in pywt.families():
         print "%s family:" % family, ', '.join(pywt.wavelist(family))
      sys.exit(1)

   if pywt.Wavelet(wtype).orthogonal == False:
      print "Warning \'%s\' is not orthogonal"%wtype

   image = recompose(prefix, levels, suffix, wtype)

   piio.write(oname, image);

if __name__ == '__main__':
    main()
