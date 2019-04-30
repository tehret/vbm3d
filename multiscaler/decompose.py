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


def decompose(image, prefix, levels, suffix, wtype, cur=0):

   print prefix+str(cur)+suffix
   piio.write(prefix+str(cur)+suffix, image)

   print image.shape
   if levels==1:
      return

   ret, _ = pywt.dwt2(image,wtype,axes=(0,1))

   decompose(ret/2.0,prefix,levels-1,suffix, wtype, cur+1)


def main():
   # verify input
   wtype = pick_option(sys.argv, 't', 'db7')
   if len(sys.argv) > 4:
      image  = piio.read(sys.argv[1])
      prefix = sys.argv[2]
      levels = int(sys.argv[3])
      suffix = sys.argv[4]
   else:
      print("Incorrect syntax, use:")
      print("  > " + sys.argv[0] + " input prefix levels suffix [-t type]")
      print("        multiscale decomposition using wavelets")
      print("        available types:")
      for family in pywt.families():
         print "%s family:" % family, ', '.join(pywt.wavelist(family))
      sys.exit(1)

   if pywt.Wavelet(wtype).orthogonal == False:
      print "Warning \'%s\' is not orthogonal"%wtype

   decompose(image, prefix, levels, suffix, wtype)


if __name__ == '__main__':
    main()
