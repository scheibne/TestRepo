
# coding: utf-8

# #This is a cleaned up version of TNO Individual Web Page.ipynb. Now showing expected exposures. 
# 
# Lynus,
# 
# To make this run, here is what you need to do.
# 
# 1) Make sure you can import all the libraries/methods in the first cell. You need to have pyOrbfit (Prof. Gerdes just / sent out instructions for downloading and installing this. You also need ccdBounds, which is on the dwgerdes GitHub / repository https://github.com/dwgerdes/tnofind
# 
# 2) Make sure you have installed ds9: http://ds9.si.edu/site/Download.html. Go with the X11 version for Mac. Test it out by going to your terminal and typing ds9. You may have to make a binary executable file.
# 
# 3) Have a directory with the following things in it: 
#     
#     1) good_2.csv
#     
#     2) exposures.csv
#     
#     3) style_content.css
# 
# 4) Go to sendObsRequest and sendSearchRequest and replace my username and password with your username and password
# 
# 5) Run this line of code: rawFileObs, rawFileSearch, observed, expected = getImageTar("good_2.csv")
# 
# 6) Go to http://desdev3.cosmology.illinois.edu:8000 (you may have enter your user name and password). Go and download the tar files and put them in the directory where everything else is stored. It may take a minute or two for the job to complete. 
# 
# 7) Run the following line of code (putting in the appropriate file path): makeIndividualWebpage('/Users/ColinS/Documents/TNOSearch','good_2.csv',expecteds, rawFileObs, rawFileSearch) 
# 
# Hopefully this works. 
# 
# -Colin
# 
# ps. You're totally right, that second break statement doesn't do anything. 
# 

# In[2]:

from __future__ import division
import ccdBounds
from pyOrbfit.Orbit import Orbit
import gzip
import glob # Lists files into a directory
import tarfile
import json
import pandas
import pylab
import ephem
import os
import numpy as np
import requests
import time
from pandas import *
from pylab import *
from astropy.io import fits
from astropy.wcs import WCS


# In[3]:

def drawObsCircle(tempfits, imgfile, regfile, ra, dec):
    with open(regfile, 'w') as fout:
        #fout.write('fk5; circle '+str(ra)+' '+str(dec)+' 6" #dash=1') This line works
        w=WCS(tempfits)
        lon, lat = np.degrees(ephem.hours(ra)),np.degrees(ephem.degrees(dec))
        pixx,pixy= w.wcs_world2pix(lon,lat,3)
        fout.write('physical; circle '+str(pixx)+' '+str(pixy)+' 6" #dash=1')
    hdu=fits.getdata(tempfits)
    h=hdu.shape[0]
    w=hdu.shape[1]
    os.system("ds9 "+tempfits+" -scale mode zscale -colorbar no -height "+str(h)+" -width "+str(w)+" -zoom to fit -region"+" "+regfile+" -saveimage png "+imgfile+" -exit")  


# In[4]:

def drawSearchEllipse(tempfits, imgfile, regfile, ra, dec, PA, a, b):
    with open(regfile, 'w') as fout:
        if a > 2:
            fout.write('fk5; ellipse '+str(ra)+' '+str(dec)+' '+str(a)+'" '+str(b)+'" '+str(PA-90)+' # dash=1')
        else:
            fout.write('fk5; box '+str(ra)+' '+str(dec)+' '+str(6)+'" '+str(6)+'" 0 #color=red dash=1')
    hdu=fits.getdata(tempfits)
    h=hdu.shape[0]
    w=hdu.shape[1]
    os.system("ds9 "+tempfits+" -scale mode zscale -colorbar no -height "+str(h)+" -width "+str(w)+" -zoom to fit -region"+" "+regfile+" -saveimage png "+imgfile+" -exit")
    


# In[5]:

def fit_orbit(df_obs):
    df_obs = df_obs.ix[['#' not in row['date'] for ind, row in df_obs.iterrows()]]   # filter comment lines
    nobs = len(df_obs)
    ralist = [ephem.hours(r) for r in df_obs['ra'].values]
    declist = [ephem.degrees(r) for r in df_obs['dec'].values]
    datelist = [ephem.date(d) for d in df_obs['date'].values]
    obscode = np.ones(nobs, dtype=int)*807
    orbit = Orbit(dates=datelist, ra=ralist, dec=declist, obscode=obscode, err=0.15)
    return orbit


# In[6]:

df_obs = read_csv('good_2.csv')
df_obs = df_obs.ix[['#' not in row['date'] for ind, row in df_obs.iterrows()]]
nobs = len(df_obs)
ralist = [ephem.hours(r) for r in df_obs['ra'].values]
declist = [ephem.degrees(r) for r in df_obs['dec'].values]
datelist = [ephem.date(d) for d in df_obs['date'].values]
obscode = np.ones(nobs, dtype=int)*807
#orbit = Orbit(dates=datelist, ra=ralist, dec=declist, obscode=obscode, err=0.15)


# In[7]:

#orb = fit_orbit(df)


# In[8]:

def compute_chip(rockra, rockdec, expra, expdec):
    '''
    Given the ra and dec of a point and of the center
    of an exposure, find the CCD containing that point.
    
    Returns a pair of the CCD name and number.
    '''
    deltara = 180/np.pi*ephem.degrees(rockra-expra).znorm  # compute difference in degrees (normalized between -180, +180)
    deltadec = 180/np.pi*ephem.degrees(rockdec-expdec).znorm  # the 180/pi is because ephem.Angle objects are natively in radians
    ccdname = 'None'
    for k in ccdBounds:
        if deltara > ccdBounds[k][0] and deltara < ccdBounds[k][1] and deltadec > ccdBounds[k][2] and deltadec < ccdBounds[k][3]:
            ccdname = k
    return ccdname, ccdNum[ccdname]



# In[9]:

#Unzips tar.gz files
def unzip_tar(tarname):
    fname = str(tarname)
    if (fname.endswith("tar.gz")):
        tar = tarfile.open(fname, 'r:gz')
        tar.extractall() ##Deleted the slash in front of c
        tar.close()
    elif (fname.endswith("tar")):
        tar = tarfile.open(fname, 'r:')
        tar.extractall()
        tar.close()
    raw_fname = fname[:36]
    return raw_fname


# In[10]:

def makeImgArray(obs_properties, flist, directory, raw_fname, exp_values):
    
    img_array=pandas.DataFrame.from_items([('expnum',[]),('refccd',[]),('refpng',[]), ('refdate',[]),('reftef',[]),('compimages',[]),('compexp',[]),('compccd',[]),('compdate',[]),('comptef',[])])
    
    def getkey(item):
            return item[0]
        
    for i in range(0,len(obs_properties)):
        
        refnum=obs_properties['expnum'][i]
        #refccd=obs_properties['ccd'][i]
        refccd=obs_properties['ccd'][i]
        refband=obs_properties['band'][i]
        refdate=obs_properties['date'][i]
        ref_ra=obs_properties['ra'][i]
        ref_dec=obs_properties['dec'][i]
        CompImgs=[]
        expnums=[]
        refpng=''
        reftef=round(float(exp_values[exp_values['expnum']==refnum]['t_eff']),3)
        
        for thumb in flist:
            
            os.chdir(directory+'/'+raw_fname+'/'+thumb)
            fit_list = glob.glob('*fits')
        
            for f in fit_list:
                h=fits.open(f)
                if h[0].header["EXPNUM"]==refnum and h[0].header['CCDNUM']==refccd:
                    reffit=f
                    refpng=thumb+'/'+reffit[:-5]+'Circ.png'
                    drawObsCircle(reffit,reffit[:-5]+'Circ.png','temp.reg',ref_ra,ref_dec)
                    for g in fit_list:
                        j=fits.open(g)
                        jnum=j[0].header['EXPNUM']
                        if (j[0].header['BAND']==refband and jnum != refnum and not (jnum in expnums)):
                            try:
                                teff=float(exp_values[exp_values['expnum']==jnum]['t_eff'])
                            except:
                                teff=0
                                print jnum
                            pnglabel=g[:-5]+'Circ.png'
                            drawObsCircle(g, pnglabel, 'temp.reg',ref_ra, ref_dec)
                            date=str(j[0].header['DATE-OBS'])
                            nicedate=date[0:4]+'/'+date[5:7]+'/'+date[8:10]+' '+date[11:19]
                            CompImgs.append([teff,thumb+'/'+pnglabel, jnum, j[0].header['CCDNUM'], nicedate])
                            expnums+=[jnum]
                        j.close()
                    flist.remove(thumb)
                    break 
                h.close()
                

        
        sortedcomps=sorted(CompImgs, key=getkey, reverse=True)
        
        tempframe=pandas.DataFrame.from_items([('expnum', refnum),('refccd',refccd), ('refpng',[ refpng]), ('refdate', [refdate]),('reftef',[reftef]),('compimages', [[x[1] for x in sortedcomps]]),('compexp',[[x[2] for x in sortedcomps]]), ('compccd',[[x[3] for x in sortedcomps]]),('compdate',[[x[4] for x in sortedcomps]]),('comptef',[[round(x[0],3) for x in sortedcomps]])])                
        img_array=img_array.append(tempframe, ignore_index=True)
        
    return img_array


# In[11]:

def makeSearchArray(obs_properties, flist, directory, raw_fname, exp_values):
    
    search_array=pandas.DataFrame.from_items([('expnum',[]),('refccd',[]),('refpng',[]), ('refdate',[]),('reftef',[]),('compimages',[]),('compexp',[]),('compccd',[]),('compdate',[]),('comptef',[])])
    
    def getkey(item):
            return item[0]
        
    for i in obs_properties.index:
        print i
        
        refnum=obs_properties['expnum'][i]
        refccd=obs_properties['ccd'][i]
        #refccd=obs_properties['ccdnum'][i]
        refband=obs_properties['band'][i]
        refdate=str(ephem.date(obs_properties['date'][i]))
        ref_ra=ephem.hours(obs_properties['can_ra'][i])
        ref_dec=ephem.degrees(obs_properties['can_dec'][i])
        ref_PA=obs_properties['PA'][i]
        ref_a=obs_properties['a'][i]
        ref_b=obs_properties['b'][i]
        CompImgs=[]
        expnums=[]
        refpng=''
        reftef=round(float(exp_values[exp_values['expnum']==refnum]['t_eff']),3)
        
        for thumb in flist:
            
            os.chdir(directory+'/'+raw_fname+'/'+thumb)
            fit_list = glob.glob('*fits')
        
            for f in fit_list:
                h=fits.open(f)
                if h[0].header["EXPNUM"]==refnum and h[0].header['CCDNUM']==refccd:
                    reffit=f
                    refpng=thumb+'/'+reffit[:-5]+'Circ.png'
                    drawSearchEllipse(reffit,reffit[:-5]+'Circ.png','temp.reg',ref_ra,ref_dec,ref_PA, ref_a,ref_b)
                    for g in fit_list:
                        j=fits.open(g)
                        jnum=j[0].header['EXPNUM']
                        if (j[0].header['BAND']==refband and jnum != refnum and not (jnum in expnums)):
                            try:
                                teff=float(exp_values[exp_values['expnum']==jnum]['t_eff'])
                            except:
                                teff=0
                                print jnum
                            pnglabel=g[:-5]+'Circ.png'
                            drawSearchEllipse(g, pnglabel, 'temp.reg',ref_ra, ref_dec, ref_PA, ref_a, ref_b)
                            date=str(j[0].header['DATE-OBS'])
                            nicedate=date[0:4]+'/'+date[5:7]+'/'+date[8:10]+' '+date[11:19]
                            CompImgs.append([teff,thumb+'/'+pnglabel, jnum, j[0].header['CCDNUM'], nicedate])
                            expnums+=[jnum]
                        j.close()
                    flist.remove(thumb)
                    break
                    break 
                h.close()
        
        sortedcomps=sorted(CompImgs, key=getkey, reverse=True)
        
        tempframe=pandas.DataFrame.from_items([('expnum', refnum),('refccd',refccd), ('refpng',[ refpng]), ('refdate', [refdate]),('reftef',[reftef]),('compimages', [[x[1] for x in sortedcomps]]),('compexp',[[x[2] for x in sortedcomps]]), ('compccd',[[x[3] for x in sortedcomps]]),('compdate',[[x[4] for x in sortedcomps]]),('comptef',[[round(x[0],3) for x in sortedcomps]])])                
        search_array=search_array.append(tempframe, ignore_index=True)
        
    return search_array


# In[12]:

def writeStatus():
    statstr=r'''<div>
                <br><br>
                <form action="action_page.php" >
                    <select name="Status">
                        <option value="Unchecked" >Unchecked</option>
                        <option value="Confirmed">Confirmed</option>
                        <option value="Rejected">Rejected</option>
                        <option value="Under Inspection">Under Inspection</option>
                    </select>
                    <input type="submit">
                        </form>
            </div>'''
    return statstr


# In[13]:

def makeHead(obs_properties, objid):
    numobs=len(obs_properties)
    headstr=r"""
    <!DOCTYPE html>
    <html>
    <head>
    <title>Candidate: """+objid+ r"""</title>
    <link href ="style_content.css" type="text/css" rel="stylesheet">
    </head>
    <body>
        <div class = "main">
            <div class = "header">
                <h1> Candidate: """+objid+r"""</h1>
                <p><a href = "homepage2.html">Home</a></p>
            </div>"""+writeStatus()+r"""
            <div class = "tab" align = "center">
                
            <table width = "700" border="1" align="center">
            <tr>
            <td align="center" colspan="10"><b>Properties</b>
            </td>
            
            <tr align = "center">
                <td>Date</td>
                <td>Ra</td>
                <td>Dec</td>
                <td>Expnum</td>
                <td>Exptime</td>
                <td>Band</td>
                <td>Ccd</td>
                <td>Mag</td>
                <td>ml_score</td>
                <td>Object ID</td>
            </tr>"""

    for i in range(0,numobs):
        headstr+= r"""        <tr align = "center">
                <td>"""+obs_properties['date'][i]+r"""</td>
                <td>"""+obs_properties['ra'][i]+r"""</td>
                <td>"""+obs_properties['dec'][i]+r"""</td>
                <td>"""+str(obs_properties['expnum'][i])+r"""</td>
                <td>"""+str(obs_properties['exptime'][i])+r"""</td>
                <td>"""+obs_properties['band'][i]+r"""</td>
                <td>"""+str(obs_properties['ccd'][i])+r"""</td>
                <td>"""+str(obs_properties['mag'][i])+r"""</td>
                <td>"""+str(obs_properties['ml_score'][i])+r"""</td>
                <td>"""+str(obs_properties['objid'][i])+r"""</td>
            </tr>"""

    headstr+="""
        </table>
            </div>"""
    return headstr



# In[14]:

def makeSearchTable(obs_properties, objid):
    headstr=r"""
            <div class = "tab" align = "center">
                
            <table width = "700" border="1" align="center">
            <tr>
            <td align="center" colspan="10"><b>Possible Exposures</b>
            </td>
            
            <tr align = "center">
                <td>Date</td>
                <td>Ra</td>
                <td>Dec</td>
                <td>Expnum</td>
                <td>Band</td>
                <td>Ccd</td>
                <td>PA </td>
                <td>a</td>
                <td>b</td>
            </tr>"""

    for i in obs_properties.index:
        print i
        headstr+= r"""        <tr align = "center">
                <td>"""+str(ephem.date(obs_properties['date'][i]))+r"""</td>
                <td>"""+str(ephem.hours(obs_properties['can_ra'][i]))+r"""</td>
                <td>"""+str(ephem.degrees(obs_properties['can_dec'][i]))+r"""</td>
                <td>"""+str(int(obs_properties['expnum'][i]))+r"""</td>
                <td>"""+str(obs_properties['band'][i])+r"""</td>
                <td>"""+str(int(obs_properties['ccd'][i]))+r"""</td>
                <td>"""+str(round(obs_properties['PA'][i],2))+r"""</td>
                <td>"""+str(round(obs_properties['a'][i],2))+r"""</td>
                <td>"""+str(round(obs_properties['b'][i],2))+r"""</td>
            </tr>"""

    headstr+=r"""
        </table>
            </div>"""
    return headstr
    


# In[15]:

def tableHeader(objid):
    tableHead = """		<div class = "img">
			<table ID="t02" >
				<tr>
					<td align="center" colspan="5"><b>Images of """+objid+"""</b>
					<p> Click on an image to expand it. </p>
					</td>
                </tr>
                    <th>Observations</th> <th>Comparison images ordered by decreasing t_eff</th> """
    return tableHead


# In[16]:

def searchHeader(objid):
    tableHead = """		<div class = "img">
			<table ID="t02" >
				<tr>
					<td align="center" colspan="5"><b>Exposures coinciding with best fit orbit of """+objid+"""</b>
					<p> Click on an image to expand it. </p>
					</td>
                </tr>
                    <th>Expected Observations</th> <th>Comparison images ordered by decreasing t_eff</th> """
    return tableHead


# In[17]:

def writeImgTable(ImgArray, directory, raw_fname):
    tableBody = ''
    for i in range(0,len(ImgArray)):
        refpng=ImgArray['refpng'][i]
        expnum=ImgArray['expnum'][i]
        refccd=ImgArray['refccd'][i]
        refdate=ImgArray['refdate'][i]
        reftef=ImgArray['reftef'][i]
        compimages=ImgArray['compimages'][i]
        compexp=ImgArray['compexp'][i]
        compccd=ImgArray['compccd'][i]
        compdate=ImgArray['compdate'][i]
        comptef=ImgArray['comptef'][i]
        tableBody+= '<tr align ="center">'
        tableBody+= '<tr><td><p></p></td></tr>'
        tableBody+= r'''<td> <u>t_eff = '''+str(reftef)+r'''</u></td>'''
        for j in range(0,len(compimages)):
            tableBody+= r'''<td><u> t_eff = '''+str(comptef[j])+'''</u></td>'''
        tableBody += '</tr>		<tr align = "center" valign="top">'
        tableBody += '			'+r''' <td><a href = "'''+directory+'/'+ raw_fname+'/'+refpng+r'''" ><img src= "'''+directory+'/'+ raw_fname+'/'+refpng+r'''" align='center' alt = "No Reference Image Available"
                width ="200" height="200" align="center" border = "1"></img> </a> </td>'''
        for j in range(0,len(compimages)):
            tableBody += '			'+r'''<td> <a href = "'''+directory+'/'+raw_fname+'/'+compimages[j]+'''"><img src="'''+directory+'/'+raw_fname+'/'+compimages[j]+'''" align='center' alt = "Example 1 pic"
                width ="200" height="200" align="center" border = "1"></img></a> </td>'''
        tableBody+= r'''</tr>     <tr align = "center">'''
        tableBody+= r'''<td>Expnum = '''+str(int(expnum))+r''', ccd = '''+ str(int(refccd))+r'''</td>'''
        for j in range(0, len(compimages)):
            tableBody+='''<td>Expnum = '''+str(compexp[j])+r''', ccd = '''+ str(compccd[j])+r'''</td>'''
        tableBody+=r'''</tr> <tr align = "center" border-bottom="1px">'''
        tableBody+=r'''<td> Date: '''+str(refdate)+r'''</td>'''
        for j in range(0, len(compimages)):
            tableBody+=r'''<td> Date: '''+str(compdate[j])+r'''</td>'''
        tableBody += "</tr> \n"
    
    return tableBody


# In[18]:

def writeComments():
    comstr = r'''<div> 
    <br><br>
    
    
    <table id="t01" width = "700" align="center" cellspacing="6" cellpadding = "2" >
        <tr>
            <th align="center" colspan="5"><b>Comments</b> </th>
        </tr>
        <tr>
            <td> 12/18/2016 </td>
        </tr>
        <tr>
        <td > TNOs are so cool! </td>
        </tr>
        <tr>
            <td > 12/15/2016 </td>
        </tr>
        <tr>
            <td> This candidate is registered with the Minor Planet Center </td>
        </tr>
    
    
    
    </table>

    </div>'''
    
    return comstr


# In[19]:

def writeLeavCom():
    formstr=r'''<div>
    <form action="gold_27.html" method="post">
        <p>Comment:</p>
        <textarea type="text" name="name"> </textarea> <br>
    
        <input type="submit">
    </form>
    </div>'''
    
    
    return formstr


# In[20]:

def buildpage(obs_properties, searchFrame, flistObs, flistSearch, directory, raw_fnameObs,raw_fnameSearch, exp_values, objid):
    
    ImgArray=makeImgArray(obs_properties, flistObs, directory, raw_fnameObs, exp_values)
    SearchArray=makeSearchArray(searchFrame,flistSearch,directory,raw_fnameSearch, exp_values)
    indivpage = makeHead(obs_properties, objid) + tableHeader(objid) + '''</tr> \n
    '''+ writeImgTable(ImgArray, directory, raw_fnameObs)+'''
        		
		
		
	</table>
		
	</div>'''+writeComments()+writeLeavCom()+ makeSearchTable(searchFrame, objid)+searchHeader(objid)+'''</tr> \n
    '''+writeImgTable(SearchArray,directory, raw_fnameSearch)+r'''
    
    </table>
    </div>'''+writeComments()+writeLeavCom()+r'''
    
    </body>

    </html>'''
    return indivpage


# In[21]:

def findExposures(candidate):
    df = read_csv(candidate)
    orb = fit_orbit(df)
    matches=DataFrame()
    for i in all_exps[all_exps['t_eff']>.3].index:
        e_ra, e_dec=all_exps['ra'][i],all_exps['dec'][i]
        pos=orb.predict_pos(all_exps['date'][i])
        pos_ra, pos_dec=pos['ra'], pos['dec']
        ccdname, ccdnum = compute_chip(pos_ra, pos_dec, e_ra, e_dec)
        if ccdnum>0:
            temp=all_exps.loc[i,['expnum','date','nite','band','t_eff']]
            temp['ccd']=ccdnum
            temp['can_ra']=pos_ra
            temp['can_dec']=pos_dec
            temp['PA']=pos['err']['PA']
            temp['a']=pos['err']['a']
            temp['b']=pos['err']['b']
            temp['isinobs']=all_exps['expnum'][i] in list(df['expnum'])
            matches=matches.append(temp,ignore_index=True)
    return matches


# In[22]:

def sendObsRequest(matches):
    ra=list(np.degrees(matches[matches['isinobs']==True]['can_ra']))
    dec=list(np.degrees(matches[matches['isinobs']==True]['can_dec']))
    bands='[g,r,i,z]'
    req='http://desdev3.cosmology.illinois.edu:8000/api?username=lzullo&password=lzu70chips&ra=%s&dec=%s&bands=%s' % (ra,dec,bands)
    submit = requests.get(req)
    return submit.json()['job']


# In[23]:

def sendSearchRequest(matches):
    ra=list(np.degrees(matches[matches['isinobs']==False]['can_ra']))
    dec=list(np.degrees(matches[matches['isinobs']==False]['can_dec']))
    bands='[g,r,i,z]'
    req='http://desdev3.cosmology.illinois.edu:8000/api?username=lzullo&password=lzu70chips&ra=%s&dec=%s&bands=%s' % (ra,dec,bands)
    submit = requests.get(req)
    return submit.json()['job']


# In[ ]:

def getImageTar(candidate):
    all_exps = read_csv('exposures.csv')
    matches=findExposures(candidate)
    rawFileObs=sendObsRequest(matches)+".tar"
    rawFileSearch=sendSearchRequest(matches)+".tar"
    observed=matches[matches['isinobs']==True]
    expected=matches[matches['isinobs']==False]
    return rawFileObs,rawFileSearch,observed,expected


# def makeIndividualWebpage(directory, obfile, tarname):
#     os.chdir(directory)
#     #obs_props=read_csv('Candidates/'+obfile)
#     obs_props=read_csv(obfile)
#     objid=obfile[:-4]
#     exposures = read_csv('exposures.csv') 
#     exp_values = exposures.ix[:,['expnum','t_eff']]
#     raw_fname = unzip_tar(tarname)
#     os.chdir(directory+'/'+raw_fname)
#     flist = glob.glob('thumbs*')    
#     page = buildpage(obs_props, flist, directory, raw_fname, exp_values, objid)
#     os.chdir(directory)
#     with open(objid+'.html', 'w') as fout:
#         fout.write(page)
#     return 'done'
# #print page

# In[ ]:

def main(directory, obfile, searchFrame, tarnameObs, tarnameSearch):
    os.chdir(directory)
    #obs_props=read_csv('Candidates/'+obfile)
    obs_props=read_csv(obfile)
    objid=obfile[:-4]
    exposures = read_csv('exposures.csv') 
    exp_values = exposures.ix[:,['expnum','t_eff']]
    raw_fnameObs = unzip_tar(tarnameObs)
    raw_fnameSearch= unzip_tar(tarnameSearch)
    os.chdir(directory+'/'+raw_fnameObs)
    flistObs = glob.glob('thumbs*')
    os.chdir(directory+'/'+raw_fnameSearch)
    flistSearch = glob.glob('thumbs*')
    page = buildpage(obs_props, searchFrame, flistObs, flistSearch, directory, raw_fnameObs, raw_fnameSearch, exp_values, objid)
    os.chdir(directory)
    with open(objid+'.html', 'w') as fout:
        fout.write(page)
    return 'done'
#print page


# In[ ]:

if __name__ == '__main__':
    
    #csv = raw_input('Enter the csv file here ')
    all_exps=read_csv('exposures.csv')
    mymatches=findExposures('good_2.csv')
    rawFileObs, rawFileSearch, observed, expecteds = getImageTar('good_2.csv') 
    main('/Users/lynuszullo/pyOrbfit','good_2.csv',expecteds, rawFileObs, rawFileSearch)
    


# In[ ]:




# In[ ]:




# makeIndividualWebpage('/Users/ColinS/Documents/TNOSearch','gold_27.csv','5ad98083-9dcb-41e3-aed2-b9bd25237430.tar.gz')

# In[ ]:




# makeIndividualWebpage('/Users/ColinS/Documents/TNOSearch','QR441.csv','8e1523a8-2194-4464-a80b-1306feadacad.tar')

# Here is an example of an input:
# makeIndividualWebpage('/Users/ColinS/Documents/TNOSearch','Fakegold.csv','5ad98083-9dcb-41e3-aed2-b9bd25237430.tar.gz')
# 

# In[ ]:

all_exps=read_csv('exposures.csv')
mymatches=findExposures('good_2.csv')


# In[ ]:

goodObs=read_csv('good_2.csv')


# In[ ]:

goodObs[goodObs['expnum']==386739]


# In[ ]:

mymatches.head()


# In[ ]:

expecteds=mymatches[mymatches['isinobs']==False]


# In[ ]:

expecteds


# In[ ]:

unzip_tar('a84925eb-f279-437c-8619-9c9b8edc2d74.tar')


# In[ ]:

os.chdir('/Users/lynuszullo/pyOrbfit'+'a84925eb-f279-437c-8619-9c9b8edc2d74')


# In[ ]:

myflist=glob.glob('thumbs*')


# In[ ]:

tempSearchArray=makeSearchArray(expecteds,myflist,'/Users/lynuszullo/pyOrbfit','a84925eb-f279-437c-8619-9c9b8edc2d74',all_exps)


# In[ ]:

tempSearchArray.head()


# In[ ]:

str(ephem.hours(expecteds['can_ra'][1]))


# In[ ]:

np.degrees(expecteds['can_ra'][1])


# In[ ]:

len(myflist)


# Potential Problem expousure numbers:
#     459985
# 230084
# 231480
# 238906
# 240458
# 245899
# 247890
# 251065
# 255875
# 261942
# 267565
# 275246
# 277592
# 354905
# 359571
# 360540
# 369023
# 372919
# 379269
# 381508
# 382518
# 389436
# 391629
# 395515
# 398240
# 400779
# 403386
# 459983
# 464794
# 466292
# 474311
# 475839
# 478352
# 483388
# 484472
# 506645
# 459985
# 459985
# 345372
# 345373
# 226647
# 228716
# 230090
# 231474
# 237666
# 242388
# 243829
# 245905
# 255881
# 258477
# 275252
# 277598
# 280306
# 345371
# 348369
# 352863
# 359564
# 367108
# 371612
# 376673
# 379263
# 381225
# 381502
# 389430
# 401525
# 475492
# 482091
# 485807
# 492427
# 494279
# 497332
# 500457
# 501985
# 506784
# 508804
# 345372
# 345373
# 345372
# 345373
# 345372
# 459985

# In[ ]:

mySearchTable=makeSearchTable(expecteds, "good_2")


# In[ ]:

str(expecteds['b'][0])


# In[ ]:

str(round(expecteds['b'][0],3))


# In[ ]:

str(ephem.date(4444444.5))


# In[ ]:

all_exps[all_exps['t_eff']>.3].index


# In[ ]:



