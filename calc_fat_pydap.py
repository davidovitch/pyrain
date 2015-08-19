# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 18:17:23 2014

@author: dave
"""

from fatigue_pydap import rainflow_astm_wrapper


rainflow_func = rainflow_astm_wrapper

sig_rf = rainflow_func(sig)
ax.hist(sig_rf, no_bins)
# place a text box in upper left in axes coords
eq = fatigue_pydap.eq_load(sig, m=m, neq=neq, rainflow_func=rainflow_func)
fat_txt = "%s\nEQ number: %d\n" % (att.name, neq)
fat_txt += "\n".join(["EQ(m=%2d):  %.4f" % (_m, eq) for _m, eq in zip(m, eq)])

QtGui.QApplication.processEvents()  # to make self.size update
ax.text(1 - (10. / self.width()), 1 - (10. / self.height()), fat_txt, transform=self.axes.transAxes, fontsize=self.font_size, horizontalalignment='right', verticalalignment='top', bbox=dict(facecolor='white'), multialignment="left")
ax.set_xlabel = ''