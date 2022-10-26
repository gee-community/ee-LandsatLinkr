# ee-LandsatLinkr

The aim of `landsatlinkr` is to make it easy to run the LandTrendr algorithm with Landsat MSS, TM, ETM, and OLI data together in Earth Engine. It assembles image collections across the five satellites that carried the MSS sensor, filters those images for quality, calculates TOA reflectance, and calculates the MSScvm cloud mask. It also allows the user to manually exclude MSS images with scan line issues or other noise. These MSS images are converted to pseudo-TM by building a relationship between MSS and TM coincident images. These converted MSS to TM images are then grouped with TM and OLI images and run through LandTrendr to track changes over time. 

## Notes

The most recent version of ee-LandsatLinkr is written for the Earth Engine Python API and executed using a Colab notebook.

As of 2022-10-25 the code is still changing frequently, so check back for updates.

Right now the best place to start is with a presentation and notebook by Annie Taylor and Justin Braaten given during a presentation at [Geo for Good '22](https://earthoutreachonair.withgoogle.com/events/geoforgood22#)

- [Presentation](https://earthoutreachonair.withgoogle.com/events/geoforgood22/watch?talk=day3-track1-talk1)
- [Slides](https://docs.google.com/presentation/d/1lBL2omCxVDRqVqu7wRoVYAHC-IgAXL-7tIHJB6ZEX6k/edit#slide=id.g941972ef10_0_0)
- [Worksheet](https://docs.google.com/document/d/1hhJhamfjfrhlPkQbd2WrYgpliilzeq9nnQj2MRPndCw/edit)
- [Notebook](https://github.com/gee-community/ee-LandsatLinkr/blob/main/landsatlinkr_workshop_template.ipynb)
