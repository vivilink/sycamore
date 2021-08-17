#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:26:11 2021

@author: linkv
"""

import Mantel
import numpy as np


def mantel(X, Y, perms = 1000):
  p = Mantel.test(X, Y, perms = perms)[1]
  if p > 3/perms: 
    return p
  else: 
    return mantel(X, Y, perms = perms * 10)


def TMRCA(tree, N):
  """
  see tskit to understand these classes: https://tskit.dev/tskit/docs/stable/python-api.html#the-tree-class
  
  - the goal of this function is to get pairwise branch lengths between all nodes of the tree
  - we go through all internal nodes and update the cells of all their descendants at once
  - the tree.time of a node gets larger backwards in time
  - the idea is to fill the TMRCA matrix with the total height of the tree and then to 
  deduct for all descendants of one node the branches from farther back in time
  for example, the terminal nodes will have all branch lengths that connect their ancestors deducted 
  -> they are only separated by the terminal branches
  Since the height of the tree can only be known by finding the oldest node, which we don't know a priori 
  we deduct the branch lengths from zero first, using the loop over nodes to find the height, and then add the height
  """

  tmrca = np.zeros([N, N])
  height = 0
  for c in tree.nodes():
    descendants = list(tree.samples(c))
    n = len(descendants)
    if(n == 0 or n == N): continue
    t = tree.time(tree.parent(c)) - tree.time(c)
    tmrca[np.ix_(descendants, descendants)] -= t
    height = max(height, tree.time(tree.parent(c))) #time returns the time of a node
  tmrca += height
  return tmrca


def make(X):
  X_ = (X+X.T)/2
  np.fill_diagonal(X_, 0)
  return X_

def diff(y, N):
  cols = np.tile(y, (N, 1))
  rows = cols.T
  buffer = cols - rows
  return np.abs(buffer)