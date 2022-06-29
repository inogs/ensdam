import argparse

import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates netcdf inputs for ensdam
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--beforedir', '-b',
                                type = str,
                                required = True,
                                help = ''' '''
                                )
    parser.add_argument(   '--afterdir', '-a',
                                type = str,
                                required = True,
                                help = ''' '''
                                )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' Path of the meshmask file'''
                                )
    parser.add_argument(   '--starttime', '-s',
                                type = str,
                                required = True,
                                help = ''' start time in yyyymmdd format'''
                                )
    parser.add_argument(   '--endtime', '-e',
                                type = str,
                                required = True,
                                help = ''' end time in yyyymmdd format'''
                                )

    return parser.parse_args()

args = argument()

import numpy as np
from commons.utils import addsep
from commons.mask import Mask
from commons.submask import SubMask
from commons.dataextractor import DataExtractor
import Sat.SatManager as Satmodule
import netCDF4
from commons.Timelist import TimeInterval, TimeList
from basins import V2 as OGS


SATDIR="/g100_work/OGS_prod100/OPA/V9C/SSPADA/GHOSH/FILES/24/SATELLITE.WEEKLY/"

TI = TimeInterval(args.starttime,args.endtime, '%Y%m%d')
TL = TimeList.fromfilenames(TI, SATDIR, "*nc", prefix="", dateformat="%Y%m%d")

RSTBEFORE=addsep(args.beforedir)
RST_AFTER=addsep(args.afterdir)
OUTDIR = addsep(args.outdir)

TheMask=Mask(args.maskfile)
coastmask = TheMask.mask_at_level(200.0)
jpk,jpj,jpi = TheMask.shape
m = 24
modvarname="P_l"
satvarname="CHL"

nFrames = TL.nTimes

nSUB = len(OGS.P.basin_list)
jpk,jpj,jpi =TheMask.shape
dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
mask0 = TheMask.cut_at_level(0)

for sub in OGS.Pred:
    SUB[sub.name]  = SubMask(sub,maskobject=mask0).mask[0,:]
    if 'atl' in sub.name: continue
    SUB['med'] = SUB['med'] | SUB[sub.name]


def dumpfile(filename,prior, posterior,obs):

    ncOUT = netCDF4.Dataset(filename,'w')
    ncOUT.createDimension("m",m)
    ncOUT.createDimension("n",nPoints)
    ncvar = ncOUT.createVariable("prior_ensemble",'f',('m','n'))
    ncvar[:] = prior
    ncvar = ncOUT.createVariable("posterior_ensemble",'f',('m','n'))
    ncvar[:] = posterior
    ncvar = ncOUT.createVariable("observations",'f',('n',))
    ncvar[:] = obs
    
    ncOUT.close()


# for ie in range(m):
#     modfile="%sRST%03d.%s.%s.nc" %(DIR,ie,timestr,modvarname)
#     Model = DataExtractor(TheMask,filename=modfile, varname=modvarname).values[0,:]
for it in range(nFrames):
    satfile=TL.filelist[it]
    Sat   = Satmodule.readfromfile(satfile,var=satvarname)
    timestr=TL.Timelist[it].strftime("%Y%m%d-12:00:00")
    
    cloudsLand = (np.isnan(Sat)) | (Sat > 1.e19) | (Sat<0)

    for isub, sub in enumerate(OGS.P.basin_list):
        outfile="%sens.%s.%s.nc" %(OUTDIR,timestr,sub.name)
        selection = ~cloudsLand & coastmask & SUB[sub.name]

        nPoints= selection.sum()
        if nPoints==0 : continue
        print(outfile)
        
        prior_ensemble=np.zeros((m,nPoints),np.float32)
        for ie in range(m):
            modfile="%sRST%03d.%s.%s.nc" %(RSTBEFORE,ie,timestr,modvarname)
            #print(modfile)
            Model = DataExtractor(TheMask,filename=modfile, varname=modvarname, dimvar=2).values
            prior_ensemble[ie,:]=Model[selection]
        
        posterior_ensemble=np.zeros((m,nPoints),np.float32)
        for ie in range(m):
            modfile="%sRST%03d.%s.%s.nc" %(RST_AFTER,ie,timestr,modvarname)
            Model = DataExtractor(TheMask,filename=modfile, varname=modvarname,dimvar=2).values
            posterior_ensemble[ie,:]=Model[selection]
            
        dumpfile(outfile, prior_ensemble, posterior_ensemble, Sat[selection])








