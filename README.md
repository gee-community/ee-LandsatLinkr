# ee-LandsatLinkr

The aim of `landsatlinkr` is to make it easy to run the LandTrendr algorithm with Landsat MSS, ETM, and OLI data together in Earth Engine. It assembles image collections across the five satellites that carried the MSS sensor, filters those images for quality, calculates TOA reflectance, and calculates the MSScvm cloud mask. It also allows the user to manually exclude MSS images with scan line issues or other noise. These MSS images are converted to pseudo-TM by building a relationship between MSS and TM coincident images. These converted MSS to TM images are then grouped with TM and OLI images and run through LandTrendr to track changes over time. 

## Module import

Include the following line at the top of every script to import the landsatlinkr library. This line is already included in the template script that you will work from in the next step. 

```js
var llr = require('users/jstnbraaten/modules:landsatlinkr/landsatlinkr.js');
```

## Pre-Processing Steps

### Prep Step 1. Create a script folder 
Create a script folder or repository to hold your landsat linkr scripts. You will ultimately create a different script for each WRS-1 tile footprint that you work with, so it is helpful to have one folder where all of your landsat linkr scripts are stored. Yours might look like this: `users/user_name/LLR_Projects`. 

### Prep Step 2. Copy the landsat linkr template script into your folder 
Add a copy of the run_llr script to your folder (there are two ways to do this):
- Open the run_llr script at this [link](https://code.earthengine.google.com/7e4b572a6c67b535f5a4aabd9dbe3f67) and save a copy to your landsat linkr project folder

OR
- Add the [landsat linkr repository](https://code.earthengine.google.com/?accept_repo=users/jstnbraaten/modules) to your EE account and save a copy of run_llr_template.js to your landsat linkr project folder

### Script Step 1. Identity MSS WRS-1 tile ID(s)
Determine the WRS-1 tile ID for your study area. WRS-1 tile IDs refer to the image footprints of the MSS sensor. To do this, open your copy of the run script (see above). We'll be editing this run script to run the Landsat Linkr workflow for your study area. 

1. In the run script, set LLR_STEP to 1 and run the script:
```js
var LLR_STEP = 1;
```
2. In the map, zoom to your study area and select its location to reveal the WRS-1 tile IDs at that location (this may take a minute or two to load)
3. Copy the WRS-1 tile ID from the pop-up window. If your study area overlaps with two or more tile footprints, you should note all of the WRS-1 tile IDs and process as many tiles that intersect the study region one at a time -- the results can be composited later. For now, select one WRS-1 tile ID to begin with (you can only process one tile ID at a time).
4. Paste this WRS-1 tile ID into your run script as the WRS_1_GRANULE variable:
```js
var WRS_1_GRANULE = '046033';
```

### Prep Step 4. Create Asset Folders
Create new Asset folders to store the intermediate and final processing results of the the Landsat Linkr workflow. To make this easy, we have provided a python notebook/script that creates the appropriate folders in your Earth Engine Assets tab.

1. Open the [LandsatLinkr Asset Manager Script](https://gist.github.com/jdbcode/36f5a04329d5d85c43c0408176c51e6d) and click Open in Colab to open the script as a python notebook
2. Run the first code block to authenticate to the EE Python API (you'll be asked to open a link with your EE account and then to copy/paste the access code)
3. Set `wrs-granule-1` to your study area's WRS-1 tile ID (same as the WRS_1_GRANULE above)
4. Set `project_dir` to 'users/your_EE_username/LandTrendr'
5. Run the second, third, and fourth code blocks (Set up project dirs and Create project dirs, the rest of the code isn't needed)
6. Check your Assets tab to see that the following folders were created:
	* users/your_EE_username/LandTrendr/your_WRS_tile_ID/WRS1_to_TM and users/your_EE_username/LandTrendr/your_WRS_tile_ID/WRS1_to_WRS2


## Running Landsat Linkr
Now that you have your script and asset folders set up, we can get started with the Landsat Linkr workflow.  

### Script Step 2. MSS WRS-1 reference image
Create an MSS WRS-1 reference image for MSS WRS1 to MSS WRS2 harmonization. This reference image will be used for spectral normalization, as some of the sensors are inconsistent. This step will create a image 'ref' in your project's asset folder, and should take about 15 minutes.

1. Set the PROJ_PATH variable to the asset folder that you created.
```js
var PROJ_PATH = 'users/your_EE_username/LandTrendr';
// Must be the same path used to create the asset folder - cannot contain / at end.
```
2. Confirm that the WRS_1_GRANULE variable is set to the WRS-1 tile ID that you identified above 
```js
var WRS_1_GRANULE = '046033';
```

3. Set the CRS variable to the geographic projection that is most relevant to your study area. If you plan on calculating areas of disturbance (or other area calculations), use a projection that preserves area such as Albers equal area conic. If there are other datasets that you'll be using in combination with this Landsat analysis, consider using the CRS of those datasets. Otherwise, you can stick with the default CRS of the EE API, which is Web Mercator (EPSG:3857).
```js
var CRS = 'EPSG:3857';
```

4. Set DOY_RANGE to indicate the days of the year that you would like to include in your analysis. Here are some examples.
To include the entire year:
```js
var DOY_RANGE = [1, 365];
```
To include the months of August and September (these vary by 1 in a leap year): 
```js
var DOY_RANGE = [213, 273];
```
5. Set MAX_CLOUD to the maximum allowable percentage of cloudiness for the included images (the default is 50%). Keep in mind that a cloud mask is applied in preparing the final image collection, so some cloudiness may be acceptable during the initial filter. 
``` js
var MAX_CLOUD = 50;
```
6. Set MAX_GEOM_RMSE to the maximum allowable spatial offset in units of a pixel (the default is one half of a pixel, or 0.5).  
```js
var MAX_GEOM_RMSE = 0.5;
```
7. Set LLR_STEP to 2 and run the script.
```js
var LLR_STEP = 2;
```
8. Once it appears, run the task `MSS-reference-image` (don't change the export settings). This should take about 15 minutes. Once the task is completed and you see an image called 'ref' in your project's asset folder, you can move onto the next step. 

### Script Step 3. Preview and filter out bad MSS images
This step will print thumbnails of all of the available MSS images based on your parameters so that you can assess image quality and remove the bad images from your image collection. This step will result in an updated `EXCLUDE_IDs` variable, and can take between 10 to 30 minutes (depending on the number of images available).

1. Set LLR_STEP to 3 and run the script.
```js
var LLR_STEP = 3;
```

2. Once the image thumbnails have loaded in the console, you'll need to inspect them for issues so that you can exclude bad images from your final collection. A bad image is one with:
- lots of discolored lines
- lines that are shifted
- thin clouds and haze
- very bright or very dark colorization compared to the others

Landsat Linkr can fix (minimize/reduce to tolerable errors):
- a few discolored lines
- thick clouds

When you find a bad image, copy it's image ID (printed just above the thumbnail) and paste it as an element in the list variable `EXCLUDE_IDS`. Image IDs should be listed as strings (in quotes) and separated by a comma. Replace any existing image IDs from the example script. 

```js
var EXCLUDE_IDS = [
  'LM20460331975128AAA04',
  ... // more image IDs
  'LM10460331975227GDS03'
 ]; 
```

If there are one or two good images for a year near the target date, it is recommended that you discard any images from the same year that are only marginally good (images with more clouds, thin clouds, a few discolored lines, etc). If there are no really good images for a given year, you can include all of the marginally good ones (Landsat Linkr attempts to mask discolored lines and clouds and uses an intra-annual medoid reduction to identify the best pixel).

This step can take some time, particularly if you have included most or all of each year (see DOY_RANGE above). Once you have previewed all of the MSS images and filtered out the bad ones (by copy/pasting those image IDs to the EXCLUDE_IDs variable), you can move onto the next step. 

### Script Step 4. Prepare WRS-1 images
This step will correct each included MSS image to the reference image `ref` that you created earlier. After running this step, all of the selected MSS imagery will be harmonized within the MSS time series. This step exports annual composite images that are based on the best pixel for each year, and adds these annual composites to the WRS1 to WRS2 asset folder that you created earlier. This step takes about one hour to run. 

1. Set LLR_STEP to 4 and run the script.
```js
var LLR_STEP = 4;
```

2. This step will create 11 tasks in the EE Task tab, one for each year of analysis. Once they appear, run each task using the default settings. This should take about an hour. 

3. Once the tasks are completed and you see 11 images named for each year in the `WRS1_to_WRS2` asset folder, you can move onto the next step. You can add these asset images to the map if you would like to check them. 

### Script Step 5. Correlate the MSS to coincident TM images
This step builds the relationship between MSS and TM images by finding MSS and TM coincident images (Landsat 4 had both sensors) and using this overlap to build a robust linear regression to link the two sensors. The resulting table of correlation coefficients (exported as an asset) contains the slope, intercept, and RMSE for each image pair within your tile (and for each band). This step requires two tasks and takes about 45 minutes to run. 

1. Set LLR_STEP to 5 and run the script.
```js
var LLR_STEP = 5;
```

2. This step will create two tasks in the EE Task tab, one for to generate correlation coefficients `mss2TmCoefCol` and one to generate an offset image `medianOffset`. Once they appear, run each task using the default settings. This should take about 45 minutes.

3. Once the tasks are completed and you see two new assets in your project's asset folder (mss2TmCoefCol and mss_offset), you can move onto the next step.

### Script Step 6. Harmonize MSS images to TM and export the corrected images
This step harmonizes the WRS-1 MSS images to the TM time series using the coefficients calculated in the previous step. Note that Landsat Linkr doesn't harmonize the WRS-2 Landsat 4 and 5 images at all because there are coincident TM images which are preferred (except for 1993 when there are few TM images). This step should take a little over an hour to run. 

1. Set LLR_STEP to 6 and run the script.
```js
var LLR_STEP = 6;
```

2. This step will create 11 tasks in the EE Task tab, one for each year of analysis. Once they appear, run each task using the default settings. This should take about an hour or longer. 

3. Once the tasks are completed and you see 11 images named for each year in the `WRS1_to_TM` asset folder, you can move onto the next step. 

### Script Step 7. Inspect time series
In the assets folder, we have a harmonized MSS collection that is ready to input into the Landtrendr algorithm in EE (add link to the EE guide here). This step collects all of the Landsat imagery from later years (ETM and OLI sensors), masks out cloudy pixels and takes the annual medoid composite. 

1. Set LLR_STEP to 7 and run the script. This creates the entire Landsat image collection and saves it as a new variable `col`.
```js
var LLR_STEP = 7;
```

2. Add the image collection to the map using the example code below and run the script. 
```js
Map.addLayer(col, {}, 'Landsat Collection');
```

3. Use the inspector to view the time series at a given point. You can check the time series chart in the Inspector window, making sure there is not a step function at or around 1984.

### Script Step 8. Run Landtrendr on the entire time series

1. Set LLR_STEP to 8 and run the script. This runs the LandTrendr algorithm on the entire Landsat image collection and saves it as a new variable `lt`.
```js
var LLR_STEP = 8;
```

2. This will add the LandTrendr result to the map.

### Script Step 9. Display greatest disturbance map

1. Set LLR_STEP to 9 and run the script.  
```js
var LLR_STEP = 9;
```
2. This runs the LandTrendr algorithm on the entire Landsat image collection and saves it as a new variable `lt`. It then runs displayGreatestDisturbance() on this result to produce two layers: 
- Magnitude of greatest disturbance (change in NDVI scaled by 100)
- Year of greatest disturbance (rainbow color ramp: 1972 in purple, 2020 in red)
