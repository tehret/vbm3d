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


def recompose(image, coarse, wtype):

   print image.shape
   
   # nothing to be done: case in which coarse is the size of image
   if coarse.shape == image.shape:
      return coarse

   coeffs2 = pywt.dwt2(image,wtype,axes=(0,1))
   if coeffs2[0].shape == coarse.shape:
      ncoarse = coarse*2.0
   elif coeffs2[0].shape[0] > coarse.shape[0] and coeffs2[0].shape[1] > coarse.shape[1]:
      ncoarse = recompose(coeffs2[0],coarse*2.0,wtype)
   else:
      print "ERROR: Are you sure that \'%s\' is the right wavelet? Sizes don't match!"%wtype
      print image.shape
      print coarse.shape
      exit()
      

   #adjust size of upsampled version
   isz = coeffs2[0].shape
   ncoarse = ncoarse[0:isz[0],0:isz[1],:]

   ncoeffs2 = (ncoarse, coeffs2[1])
   ret = pywt.idwt2(ncoeffs2,wtype,axes=(0,1))

   return ret



def main():
   # verify input
   wtype = pick_option(sys.argv, 't', 'db7')
   if len(sys.argv) > 3:
      image  = piio.read(sys.argv[1])
      coarse = piio.read(sys.argv[2])
      oname  = sys.argv[3]
   else:
      print("Incorrect syntax, use:")
      print("  > " + sys.argv[0] + " image coarse result [-t type]")
      print("        2 scale recomposition using wavelets")
      print("        available types:")
      for family in pywt.families():
         print "%s family:" % family, ', '.join(pywt.wavelist(family))
      sys.exit(1)

   if pywt.Wavelet(wtype).orthogonal == False:
      print "Warning \'%s\' is not orthogonal"%wtype

   res = recompose(image, coarse, wtype)

   piio.write(oname, res);

if __name__ == '__main__':
    main()
