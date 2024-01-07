# Introduction

We built some terrain and cosmonuclide accumulation models to interpret a cosmonuclide depth profile from a vertical core taken at the top of a [tor](https://en.wikipedia.org/wiki/Tor_(rock_formation)) close to [Ben Avon](https://en.wikipedia.org/wiki/Ben_Avon).

## Core

A 3.5 m vertical core was taken inside a pothole at the top of the granitic tor. An additional "Top sample" was taken at the top surface of the tor, 123 cm above the base of the pothole.

![DSC_0156_sampling](https://github.com/angelrodes/Ben-Avon/assets/53089531/7beb17ff-9c9e-44e4-a5be-2e8e80ae10bb)\
_Core-sampling at the top of the tor._

# Methods

## Cosmogenic data

### Concentrations

<sup>10</sup>Be, <sup>26</sup>Al, and <sup>21</sup>Ne concentrations were measured in quartz from the top and core samples.

![image](https://github.com/angelrodes/Ben-Avon/assets/53089531/c6ef81bc-044b-4fdc-b6b8-21ae6858343a)\
_Sample concentrations are depicted with red stars, and one-sigma uncertainties are represented by left and right-pointing triangles._

- Both <sup>10</sup>Be and <sup>26</sup>Al profiles show an exponential distribution with an apparent attenuation length of c. 0.5 m, similar to the spallation attenuation length in granites (ρ=2.6 g/cm<sup>2</sup>), which is Λ = 160/2.6 ≈ 60 cm. This is compatible with a near-surface accumulation, where spallation is the main process that produce cosmogenic isotopes.
- The <sup>21</sup>Ne profile also show a similar "near-surface" attenuation, excluding the top sample, that show a <sup>21</sup>Ne concentration lower than the highest sample from the core.
- The <sup>21</sup>Ne profile is more scattered than the <sup>10</sup>Be and <sup>26</sup>Al profiles and do not converge to a 0 concentration at depth. This is a typical distribution of stable cosmonuclides profiles, as the quartz accumulates <sup>21</sup>Ne since its formation. This apparent "extra" <sup>21</sup>Ne concentration of c. 20 Matoms/g correspond, either to nucleogenic <sup>21</sup>Ne, either to <sup>21</sup>Ne accumulated several meters below the surface by muonic radiation.
- [<sup>26</sup>Al]/[<sup>10</sup>Be] ratios are between 6.3 and 10.4. These values are around surface production rate ratios (~7). These non-depleted ratios indicate that most of <sup>10</sup>Be and <sup>26</sup>Al atoms formed recently compared to their half-lives (1.4 and 0.7 Ma respectively).
- [<sup>21</sup>Ne]/[<sup>10</sup>Be] ratios are between one and two orders of magnitude higher than surface production rates (~4), even considering a non-cosmogenic [<sup>21</sup>Ne] of 22 Matoms/g (the lowest <sup>21</sup>Ne concentration). As the [<sup>26</sup>Al]/[<sup>10</sup>Be] ratios are not depleted and the [<sup>10</sup>Be]/[<sup>21</sup>Ne] ratios are highly depleted, we should assume that a high proportion of the <sup>21</sup>Ne atoms formed much before than most <sup>10</sup>Be and <sup>26</sup>Al atoms compared to their half-lives.

This dataset suggests that most of the <sup>21</sup>Ne formed near the surface (near-surface attenuation distribution) before most of the <sup>10</sup>Be and <sup>26</sup>Al, and then the quartz was shielded from comic radiation enough time for most of the <sup>10</sup>Be and <sup>26</sup>Al to decay. The recent near-surface exposure created most of the <sup>10</sup>Be and <sup>26</sup>Al (natural ratios) and a small fraction of the total <sup>21</sup>Ne concentrations.

### Surface production rates

We calculated the apparent <sup>10</sup>Be, <sup>26</sup>Al, and <sup>21</sup>Ne ages of the top sample using [The online calculators formerly known as the CRONUS-Earth online calculators v.3](https://hess.ess.washington.edu/math/v3/v3_age_in.html). We used the LSDn ages and the spreadsheet from the [Average cosmogenic production rate calculator](https://github.com/angelrodes/average_cosmogenic_production_rate_calculator) to calculate the surface production rates for <sup>10</sup>Be, <sup>26</sup>Al, and <sup>21</sup>Ne by spallation and muons, and their corresponding attenuation lengths, according to the approximations shown in [Rodés (2021)](https://doi.org/10.3390/geosciences11090362).

## Digital elevation model

As the samples are taken from inside an odd-shaped landform and the samples are exposed to cosmic rays that cross the tor in high variety of angles and rock thicknesses, a digital elevation model (DEM) is needed to calculate how much cosmic radiation the samples received during the geological history of the tor.

### Core location

The GPS data taken in the field during the sampling campaign locates the core at the geographical coordinate 57.109150 N, -3.462608 E. However, according to [Google Maps](https://maps.app.goo.gl/z3deD5SJwoSBNCmq9) and the digital terrain model (DMT) used here, the top of the tor is located 9.7 m to the north of this location. Therefore, we used the coordinates [ 311510 , 802980.7 ] in the British National Grid for the location of the samples.

![coordinates_correction](https://github.com/angelrodes/Ben-Avon/assets/53089531/d318bd7f-49ba-442c-b67b-8e7072e1aeb3)\
_Original GPS coordinates on Google Earth and distance to the top of the tor._

![coordinates_correction_DEM](https://github.com/angelrodes/Ben-Avon/assets/53089531/b7bad8fb-0f06-447a-a1b6-5c781e9e20db)\
_Original and corrected coordinates on the DMT used._

### Sources of elevation data

In order to simulate the cosmic irradiation of the core samples, we need a DEM of the terrain at the area where the tor is located, but also a very detailed surface model of the tor. To achieve this, we combined the following sources of terrain data:

- 5-m Digital Model of the Terrain from [Digimap](https://digimap.edina.ac.uk/).
- A 6-meter wide E-W cross-section measured in the field during the sampling campaign.
- A picture showing the S-N contour of the tor. This contour was digitalised using g3data software ([Frantz, 2000](https://manpages.ubuntu.com/manpages/jammy/man1/g3data.1.html)).

![DSC_0136_g3data_S-N](https://github.com/angelrodes/Ben-Avon/assets/53089531/0fa6dd2c-26c4-4b48-a261-534f4f6220df)\
_Digitalised picture using G3data Graph Analyser. The person at the right of the picture was used as scale._

![Elevation_data_sources](https://github.com/angelrodes/Ben-Avon/assets/53089531/e10d690e-736e-477f-9f4a-bd79f41bcc8b)\
_DMT from [Digimap](https://digimap.edina.ac.uk/) and digitalised data from field notes and one picture. Red stars depict sample positions._

### Radial DEM

Using all the topographic data, a polar mesh centred at the core position was generated. 100 azimuths are equally spaced every 3.6 degrees. Radial distances are logarithmic-like spaced to get a good resolution of the samples' shielding.

![3-D_DEM](https://github.com/angelrodes/Ben-Avon/assets/53089531/d0edf8c9-098f-4983-add8-4e2303495450)\
_Radial DEM generated from the three sources of topographic data. The left image shows the full DEM and the right image is restricted to the tor area. Sample locations, partially hidden by the grid, are depicted by red stars._

### Tor areas

Tor areas were defined by areas protruding more than 0.5 above a smoothed surface, interpolated from areas far from the centre of the grid.

![DEM_and_satellite](https://github.com/angelrodes/Ben-Avon/assets/53089531/55693274-2504-478d-9250-f2f3a9c529b7)\
_Defined tor areas compared to the satellite image from the same area._

## Erosion models

Three erosion models were build in order to test three different scenarios for the formation of the tor:

- Model A: This model considers the tor as an antecedent landform. The relief migrates down with time, with a constant erosion rate for all points of the DEM. In this model, the tor is located in an area that was originally more prominent than the surrounding areas. The surface position is calculated by adding vertical lowerings to the current DEM.
- Model B: This model simulates that the tor is "exhumed" by differential erosion from an original smooth surface. In this model, the tor formed due to its differential erosionability (e.g. harder granite or fewer joints). The surface position is calculated by adding vertical lowerings to ta smooth version of the current DEM excluding the tors. In the defined tor areas, the surface is calculated as the maximum between the smoothed DEM and the current DEM.
- Model C: This model assumes that the tor is formed by lateral erosion. In this model, the tor is conserved due to its distance to the areas where ledge migration initiated. The surface positions are calculated by selecting the elevations of the current DEM and of an elevated smoothed DEM covering the tor. The selection is based in the distance to the centre of the tor and the "horizontal lowering" for each evolution stage.

![ModelA](https://github.com/angelrodes/Ben-Avon/assets/53089531/279f01f0-ed9c-4cc6-9875-0961158865cc)\
_Model A: Antecedent landform._

![ModelB](https://github.com/angelrodes/Ben-Avon/assets/53089531/626c3d7b-61ec-4bc9-a453-c92bf708d0a0)\
_Model B: Exhumed landform._

![ModelC](https://github.com/angelrodes/Ben-Avon/assets/53089531/e477cf28-28f9-4cff-8306-32b47755a3bf)\
_Model C: Ledge migration._

## Shielding factors

In order to calculate the production of cosmogenic isotopes in samples below the surface, we need to calculate how much of the cosmic radiation is shielded in all directions by the rock between the sample and the surface. [Balco (2014)](https://doi.org/10.1016/j.quageo.2013.12.002) described a method for estimating cosmic-ray shielding by oddly shaped objects. Here we developed a slightly simpler method that require affordable computing times for large DEMs. Instead of randomizing the trajectory of the cosmic rays that hit the sample, our method calculates the attenuation of the cosmic rays that cross any point of the DEM above the sample position. As we built a radial DEM with logarithmic-like distance distribution, our method simulates the irradiation of the sample with a homogeneous distribution of zenith and azimuthal angles for the cosmic rays.

![Shielding_method](https://github.com/angelrodes/Ben-Avon/assets/53089531/36446801-ed40-46c7-aa93-4e46375ec4ce)\
_Shielding calculations._

For each trajectory, the model calculates the distance (d) between the sample and the surface and the zenith angle (φ). The attenuation of the cosmic ray for this trajectory is calculated as e<sup>-d*ρ/Λ</sup>, where ρ is the density of the granite (2.6 kg/L) and Λ is the attenuation length of the cosmic radiation (160 g/cm<sup>2</sup> for spallation and 850, 5000, and 500 g/cm<sup>2</sup> for muons).

All the attenuations calculated are weighted by N<sub>ΔΦ</sub> · _sin_(φ)<sup>2.3</sup>, where N<sub>ΔΦ</sub> is the number of cosmic ray for each of the 100 zenith angle groups (one every 0.9°), and _sin_(φ)<sup>2.3</sup> is the relative contribution of the cosmic radiation according to [Gosse & Phillips (2001)](https://doi.org/10.1016/S0277-3791(00)00171-2). The final shielding factor is calculated as 1 minus the weighted average attenuation.

In contrast to Balco's method, holes that make cosmic rays to cross the surface 3 or more times are not simulated, and their contributions are considered "in average" according to the different distances between the surface and the sample. As most of the samples are under the surface, we don't expect many close-to-vertical rays, the ones that contribute more, to cross the surface 3 or more times, so we consider that this approximation provide accurate shielding factors.

![Shielding_profile_calculations](https://github.com/angelrodes/Ben-Avon/assets/53089531/bc829d86-11c7-445f-a267-8e6aec6b89cf)\
_Shielding factors were calculated for each DEM, each sample, and each attenuation length considered for spallation and muons. These graphs show the shielding calculations for the current DEM._

![Shielding_models](https://github.com/angelrodes/Ben-Avon/assets/53089531/8e9ab89e-6e92-4629-8bd4-c66d94fcb323)\
_Distribution of shielding factors of each DEM generated for each erosion model and the different attenuation lengths._

## Glacial model

As the <sup>10</sup>Be, <sup>26</sup>Al, and <sup>21</sup>Ne concentrations indicate that this landform has been through a complex exposure history (at least exposure-burial-exposure), we tested the model used in [Sudgen _et al._ (2017)](https://doi.org/10.1016/j.epsl.2017.04.006) in this dataset. This model simulates the effect of intermittent glacial cover of the surface by considering that the samples are exposed to near-surface cosmic radiation during ice-free stages (interglacials) and shielded from cosmic radiation during ice-covered periods (glaciations). The distribution of glacial-interglacial stages is based on the δ<sup>18</sup>O curve and the δ<sup>18</sup>O threshold simulated for each model. The calculator simulates and finds both the age of the landform and the δ<sup>18</sup>O threshold that best fit the cosmogenic depth profiles. The composite δ<sup>18</sup>O curve was taken from the [NUNAIT code](https://github.com/angelrodes/NUNAIT) described in [Rodés (2021)](https://doi.org/10.3390/geosciences11090362), as it is optimised for cosmogenic accumulation calculations.

![Climate_curve](https://github.com/angelrodes/Ben-Avon/assets/53089531/bd3a4804-6bd7-403f-9eda-13d54c004971)\
_Composite climate proxy used to simulate exposure and burial stages through the glacial history of the surface._

# Results

## Erosion Models

![results_erosion_models](https://github.com/angelrodes/Ben-Avon/assets/53089531/f86528de-261f-4fd6-99fe-7b297687c949)\
_Best fit of each of the tested model (blue lines) and sample concentrations (red stars). Vertical and horizontal exhumation rates of the best fitting models are shown in bold._

- As expected, all these models fail to simulate the high <sup>21</sup>Ne concentrations.
- Models A and B (antecedent and exhumed landform) simulate the <sup>10</sup>Be and <sup>26</sup>Al depth profiles much better than model C (ledge migration).
- Both best fitting exhumation rates of models A and B suggest that the top of the tor accumulated <sup>10</sup>Be and <sup>26</sup>Al equivalent to an apparent exposure time of 210 ka.

## Glacial model

![results_glacial_model](https://github.com/angelrodes/Ben-Avon/assets/53089531/fbe35af9-1d3c-4a6e-bea6-81ab1831e886)\
_Upper graphs: best fit of the tested model (blue lines) and sample concentrations (red stars). Lower graph: The blue star show the formation of the landform and the blue line its evolution respect the δ<sup>18</sup>O curve. When the blue line is above the black δ<sup>18</sup>O curve, the landform is exposed to cosmic radiation; and when the blue line is below the black δ<sup>18</sup>O curve, the landform does not receive any cosmic radiation due to complete ice-shielding. Note that, in contrast to nunataks, higher δ<sup>18</sup>O thresholds correspond to lower elevations when simulating ice-cover by mid-latitude glaciers. This is because lower areas are less likely covered by valley glaciers and higher areas are less likely covered by polar ice-sheets in nunataks._

Although the model dos not fit the <sup>10</sup>Be and <sup>26</sup>Al profiles as good as the erosion models, this is the only one that succesfully emulates the shape of the <sup>21</sup>Ne profile while fitting the <sup>10</sup>Be profile. Assuming a possible non-cosmogenic <sup>21</sup>Ne contributon of ~20 Matoms/g, this model fits the <sup>21</sup>Ne dataset sufficiently.

The best fitting complex-exposure glacial model implies a 6.8 Ma history:

  - From 6.8 to ~4 Ma, the landform was mostly exposed to cosmic radiation near the surface.
  - Since ~4 Ma, the landform was mostly shielded from cosmic radiation due to ice-cover.

The best fitting δ<sup>18</sup>O threshold is much lower than current δ<sup>18</sup>O values. In this mid-latitude context, the best fitting δ<sup>18</sup>O threshold correlate with a much higher elevation that today should be covered by a glacier. However, glacial isostasy could explain long term subsidence of the Cairngorms, slightly moving the elevation of the tor to lower (less glaciated) conditions during the Quaternary. This vertical migration of the surface is not simulated in the glacial model. The Quaternary subsidence of the Cairgnoms could have allowed this area to be ice-free during part of the last few glacialy cycles, permitting the accumulation of the recent <sup>10</sup>Be and <sup>26</sup>Al.

# Conclusions

- Fitted erosion models suggest that the tor formed by differential erosion or by epigenesis (inherited relief) rather than lateral erosion.
- Fitted glacial model indicate that the tor was exposed near the surface for at least 2 Ma and then covered by ice during part of the Pliocene and most of the Quaternary. However, the tor is in an area that has been ice-free at least since the Yonger Dryas. Also, the <sup>10</sup>Be and <sup>26</sup>Al data and the erosion models indicate recent near-surface conditions for c. 200 ka.
- Two scenarios could explain the conditions needed to achieve long periods of shielding since the Pliocene in a place that is now ice-free and exposed to cosmic radiation:
  - The tor was on the bed of a glacial valley at some point during the Pliocene and was covered by a thick till deposit that disappeared during the Quaternary due to lateral migration of the glacial valley.
  - The tor formed at a higher elevation during the Pliocene and covered by ice through most Plio-Quaternary glaciations, while the whole area sank due to glacial isostatic adjustment.


