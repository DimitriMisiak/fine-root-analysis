#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 22:11:20 2019

@author: misiak
"""

import scipy.stats as st
import numpy as np

class double_norm(st.rv_continuous):
    """ Double Gaussian distribution. """

    def _cdf(self, x, f, loc1, scale1, loc2, scale2):
        cdf1 = (1-f) * st.norm.cdf(x, loc=loc1, scale=scale1)
        cdf2 = f * st.norm.cdf(x, loc=loc2, scale=scale2)
        cdf = cdf1 + cdf2
        return cdf

    def _pdf(self, x, f, loc1, scale1, loc2, scale2):
        pdf1 = (1-f)* st.norm.pdf(x, loc=loc1, scale=scale1)
        pdf2 = f* st.norm.pdf(x, loc=loc2, scale=scale2)
        pdf = pdf1 + pdf2
        return pdf

    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
#        for arg in args:
#            cond = np.logical_and(cond, (np.asarray(arg) > 0))
        return cond

    def _fitstart(self, data):
        mu01, mu02 = np.quantile(data, [0.10, 0.90])
        sig01 = sig02 = abs(mu01 - mu02)/50
        p0_light = [0.5, mu01, sig01, mu02, sig02]
        p0 = np.append(p0_light, [0, 1])
        return tuple(p0)

class fid_mixture(st.rv_continuous):
    """ Double Gaussian distribution plus uniform distribution """

    def _cdf(self, x, f, loc1, scale1, loc2, scale2, fu, loc3, scale3):
        cdf1 = (1-fu) * (1-f) * st.norm.cdf(x, loc=loc1, scale=scale1)
        cdf2 = (1-fu) * f * st.norm.cdf(x, loc=loc2, scale=scale2)
        cdf3 = fu * st.uniform.cdf(x, loc=loc3, scale=scale3)
        cdf = cdf1 + cdf2 + cdf3
        return cdf

    def _pdf(self, x, f, loc1, scale1, loc2, scale2, fu, loc3, scale3):
        pdf1 = (1-fu) * (1-f)* st.norm.pdf(x, loc=loc1, scale=scale1)
        pdf2 = (1-fu) * f* st.norm.pdf(x, loc=loc2, scale=scale2)
        pdf3 = fu * st.uniform.pdf(x, loc=loc3, scale=scale3)
        pdf = pdf1 + pdf2 + pdf3
        return pdf

    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
#        for arg in args:
#            cond = np.logical_and(cond, (np.asarray(arg) > 0))
        return cond
    
    def _fitstart(self, data):
        mu01, mu02 = np.quantile(data, [0.10, 0.70])
        mu03 = np.mean(data)
        sig01 = sig02 = abs(mu01 - mu02)/50
        sig03 = abs(mu01 - mu02)
        p0_light = [0.50, mu01, sig01, mu02, sig02, 0.05, mu03, sig03]
        p0 = np.append(p0_light, [0, 1])
        return tuple(p0)
    
    
class white_norm(st.rv_continuous):
    """ Double Gaussian distribution plus uniform distribution """

    def _cdf(self, x, fu, loc1, scale1, loc3, scale3):
        cdf1 = (1-fu) * st.norm.cdf(x, loc=loc1, scale=scale1)
        cdf3 = fu * st.uniform.cdf(x, loc=loc3, scale=scale3)
        cdf = cdf1 + cdf3
        return cdf

    def _pdf(self, x, fu, loc1, scale1, loc3, scale3):
        pdf1 = (1-fu) * st.norm.pdf(x, loc=loc1, scale=scale1)
        pdf3 = fu * st.uniform.pdf(x, loc=loc3, scale=scale3)
        pdf = pdf1 + pdf3
        return pdf

    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
#        for arg in args:
#            cond = np.logical_and(cond, (np.asarray(arg) > 0))
        return cond
    
    def _fitstart(self, data):
        mu01 = np.median(data)
        mu03 = np.min(data)
        sig01 = np.std(data)
        sig03 = np.max(data) - np.min(data)
        p0_light = [0.050, mu01, sig01, mu03, sig03]
        p0 = np.append(p0_light, [0,1])
        return tuple(p0)
            