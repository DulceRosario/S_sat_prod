import sys
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import geopandas as gpd
from datetime import datetime
from goes_tools import radiance_to_bt,get_like_index
from dotenv import dotenv_values,load_dotenv
import  os

load_dotenv()
inputPath=os.getenv("DOWNLOADPATH")
outPath=os.getenv("PROD_PATH")
def prod_airmass(date):

    print("******** AIRMASS  ************* ")
    prodPath="airmass/"
    dominio_1="rep"
    dominio_2="came"
    title_dom_1="Airmass RGB - México "
    title_dom_2="Airmass RGB - Centro México "
    out_file_dom_1=outPath+dominio_1+"/"+prodPath+"airmass_"+date+"_"+dominio_1+".png"
    out_file_dom_2=outPath+dominio_2+"/"+prodPath+"airmass_"+date+"_"+dominio_2+".png"

    d1_file=inputPath+"goes16.abi-"+date+"-C08_2km.tif"
    d2_file=inputPath+"goes16.abi-"+date+"-C10_2km.tif"
    d3_file=inputPath+"goes16.abi-"+date+"-C12_2km.tif"
    d4_file=inputPath+"goes16.abi-"+date+"-C13_2km.tif"

    date_time_obj = datetime.strptime(date, '%Y.%m%d.%H%M')

    D1_dataset= xarray.open_dataset(d1_file, engine="rasterio")
    D1=(D1_dataset.Rad.sel(band=1).values) 
    ##########################CANAL8########################### 
    fk1=5.03614e+04
    fk2=2.32657e+03
    bc1=2.12504
    bc2=0.99541
    D1T=radiance_to_bt(D1,fk1,fk2,bc1,bc2)-273

    D2_dataset= xarray.open_dataset(d2_file, engine="rasterio")
    D2=D2_dataset.Rad.sel(band=1).values
    ##########################CANAL10########################## 
    fk1=3.00925e+04
    fk2=1.95961e+03
    bc1=0.06984
    bc2=0.99983
    D2T=radiance_to_bt(D2,fk1,fk2,bc1,bc2)-273

    D3_dataset= xarray.open_dataset(d3_file, engine="rasterio")
    D3=D3_dataset.Rad.sel(band=1).values
    ##########################CANAL12########################## 
    fk1=1.34382e+04
    fk2=1.49784e+03
    bc1=0.10861 
    bc2=0.99966
    D3T=radiance_to_bt(D3,fk1,fk2,bc1,bc2)-273

    D4_dataset= xarray.open_dataset(d4_file, engine="rasterio")
    D4=D4_dataset.Rad.sel(band=1).values
    ##########################CANAL13########################## 
    fk1=1.07364e+04
    fk2=1.38986e+03
    bc1=0.13445
    bc2=0.99955
    D4T=radiance_to_bt(D4,fk1,fk2,bc1,bc2)-273


    # RGB Components
    RT = D1T - D2T 
    GT = D3T - D4T
    BT = D1T

    # Minimuns and Maximuns
    Rmin = -26.2
    Rmax = 0.6
    Gmin = -43.2
    Gmax = 6.7
    Bmin = -29.25
    Bmax = -64.65

    RT[RT<Rmin] = Rmin
    RT[RT>Rmax] = Rmax
    GT[GT<Gmin] = Gmin
    GT[GT>Gmax] = Gmax
    BT[BT<Bmax] = Bmax
    BT[BT>Bmin] = Bmin

    # Choose the gamma
    gamma_R = 1
    gamma_G = 1
    gamma_B = 1
    # Normalize the data
    R_bright = ((RT - Rmin) / (Rmax - Rmin)) ** (1/gamma_R)
    G_bright = ((GT - Gmin) / (Gmax - Gmin)) ** (1/gamma_G)
    B_bright = ((BT - Bmin) / (Bmax - Bmin)) ** (1/gamma_B) 

    RGB = np.dstack([R_bright, G_bright, B_bright])

    bounds = (D1_dataset.x.data[0],D1_dataset.x.data[-1],D1_dataset.y.data[-1],D1_dataset.y.data[0])
    gdf = gpd.read_file('archivos_extra/destdv1gw/destdv1gw.shp')


    fig = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.3,color='white',linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'size': 6, 'color': 'white'}
    gl.xlabel_style = {'size': 6, 'color': 'white'}
    gl.ypadding = -10 
    gl.xpadding = -10 
    ax.coastlines(color='white', linewidth=0.2)
    ax.imshow(RGB, extent = bounds, origin = 'upper')
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.08)
    plt.title(title_dom_1+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9, x=0.43, y=0.95)
    plt.savefig(out_file_dom_1,bbox_inches = 'tight', pad_inches = 0 )

    latitudes=D1_dataset.y.values
    longitudes=D1_dataset.x.values
    max_lat= get_like_index(latitudes,21.95 )
    min_lat=get_like_index(latitudes,15.3468056)
    min_lon=get_like_index(longitudes,-102.96770000)
    max_lon=get_like_index(longitudes,-93.19)

    bounds2 = (longitudes[min_lon],longitudes[max_lon], latitudes[min_lat],latitudes[max_lat])

    fig = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.imshow(RGB[max_lat:min_lat,min_lon:max_lon], extent = bounds2, origin = 'upper')
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.15)
    ax.set_xlim(longitudes[min_lon], longitudes[max_lon])
    ax.set_ylim(latitudes[min_lat], latitudes[max_lat])
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.2,color='white',linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'size': 6, 'color': 'white'}
    gl.xlabel_style = {'size': 6, 'color': 'white'}
    gl.ypadding = -10 
    gl.xpadding = -10 
    plt.title(title_dom_2+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9, x=0.43, y=0.95)

    plt.savefig(out_file_dom_2,bbox_inches = 'tight',pad_inches = 0 )