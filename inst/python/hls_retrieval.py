import requests
import json
import time
import os
from pandas.io.json import json_normalize
from collections import defaultdict
import pandas as pd
import zipfile, io

values={"aoi": '{"type":"Feature","geometry":{"type":"Polygon","coordinates":[[[-92.329005, 41.6850758], [-92.3290479, 41.681366], [-92.3240805, 41.6813179], [-92.3255289, 41.6841784], [-92.329005, 41.6850758]]]}}',
"Band": "['NDVI']",
"Enddate": "3/8/2019",
"Startdate": "3/2/2019",
"legendtype": "Relative",
"qafilter": "1",
"satellite": "Landsat,Sentinel",
"resolution":"0.00001"}
    
headers = {'Content-Type':'application/x-www-form-urlencoded'}
           
def harmonized_landsat_data(values, headers):
    try:
        url = 'https://ag-analytics.azure-api.net/hls-service/'
        response = requests.post(url, data = values, headers = headers).json()
        return response
    
    except Exception as e:
        raise e

hlsresponse = harmonized_landsat_data(values, headers)
print(hlsresponse[0].keys())
