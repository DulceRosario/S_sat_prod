import rioxarray
import xarray
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import geopandas as gpd
from datetime import datetime
from goes_tools import creapaleta, radiance_to_bt, temptoindex,get_like_index
from PIL import Image, ImageDraw, ImageFont, ImageColor, ImagePalette,ImageEnhance
import matplotlib.pyplot as plt
from dotenv import dotenv_values,load_dotenv
import  os

load_dotenv()
inputPath=os.getenv("DOWNLOADPATH")
outPath=os.getenv("PROD_PATH")

def prod_infrarojo(date):
    print("Chanel_13 :",date)
    prodPath="infrared/"
    dominio_1="rep"
    dominio_2="came"
    title_dom_1="Temperatura tope Nube - México"
    title_dom_2="Temperatura tope Nube - Centro México "
    out_file_dom_1=outPath+dominio_1+"/"+prodPath+"infrared_"+date+"_"+dominio_1+".png"
    out_file_dom_2=outPath+dominio_2+"/"+prodPath+"infrared_"+date+"_"+dominio_2+".png"
   

    d4_file="Images/goes16.abi-"+date+"-C13_2km.tif"
    D4_dataset= xarray.open_dataset(d4_file, engine="rasterio")
    D4=D4_dataset.Rad.sel(band=1).values
    bounds = (D4_dataset.x.data[0],D4_dataset.x.data[-1], D4_dataset.y.data[-1], D4_dataset.y.data[0])
    
    fk1=1.07364e+04
    fk2=1.38986e+03
    bc1=0.13445
    bc2=0.99955

    D4T=radiance_to_bt(D4,fk1,fk2,bc1,bc2)-273
    temp_min=-80
    temp_max = np.nanmax(D4T)
    temp_bar=-30

    minTempIdx = int(temptoindex(temp_bar, temp_max, temp_min)) 
    temp = ((temp_max - D4T)*256)/(temp_max - temp_min + .0001)
    p = creapaleta(minTempIdx)
    PIL_image = Image.fromarray(temp.astype('uint8'))
    enhancer = ImageEnhance.Sharpness(PIL_image)

    im_output = enhancer.enhance(1.3)
    im_output.putpalette(p)
    gdf = gpd.read_file('archivos_extra/destdv1gw/destdv1gw.shp')

    fig = plt.figure(figsize=(10.24,9.1),)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.3, color='white',linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'size': 7, 'color': 'white','weight': 'bold'}
    gl.xlabel_style = {'size': 7, 'color': 'white','weight': 'bold'}
    gl.ypadding = -10 
    gl.xpadding = -10 
    ax.coastlines(color='white', linewidth=0.2)

    ax.imshow(im_output, extent = bounds, origin = 'upper')
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.1)
    #plt.title(title_dom_1+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9,x=0.43, y=0.95)
    plt.savefig(out_file_dom_1,bbox_inches = 'tight', pad_inches = 0,dpi=250)


    latitudes=D4_dataset.y.values
    longitudes=D4_dataset.x.values
    max_lat= get_like_index(latitudes,21.95 )
    min_lat=get_like_index(latitudes,15.3468056)
    min_lon=get_like_index(longitudes,-102.96770000)
    max_lon=get_like_index(longitudes,-93.19)
    bounds2 = (longitudes[min_lon],longitudes[max_lon],latitudes[min_lat],latitudes[max_lat])

    temp2=temp[max_lat:min_lat,min_lon:max_lon]
    p = creapaleta(minTempIdx)
    PIL_image2 = Image.fromarray(temp2.astype('uint8'))
    enhancer = ImageEnhance.Sharpness(PIL_image2)
    im_output2 = enhancer.enhance(1.3)
    im_output2.putpalette(p)

    fig = plt.figure(figsize=(10.24,9.1))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.imshow(im_output2, extent = bounds2, origin = 'upper')
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.3)
    ax.set_xlim(longitudes[min_lon], longitudes[max_lon])
    ax.set_ylim(latitudes[min_lat], latitudes[max_lat])
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.35,color='white',linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'size': 7, 'color': 'white','weight': 'bold'}
    gl.xlabel_style = {'size': 7, 'color': 'white','weight': 'bold'}    
    gl.ypadding = -10 
    gl.xpadding = -10 
    #plt.title(title_dom_2+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9,x=0.43, y=0.95)

    plt.savefig(out_file_dom_2,bbox_inches = 'tight', pad_inches = 0,dpi=120)  


