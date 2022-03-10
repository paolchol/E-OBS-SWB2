# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 12:08:18 2022

@author: paolo
"""

def some_func():
    print('this is some_func')
    print(__name__)

if __name__ == '__main__':
    some_func()
    print('this is script test1')

if __name__ == 'test1':
    some_func()
    print('this is script test1')
    print('newstuff')
