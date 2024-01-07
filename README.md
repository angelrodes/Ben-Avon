# Introducction

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
_Sample concentrations are depicted with red stars and one-sigma uncertainties are represented by left and right-pointing triangles._

- Both <sup>10</sup>Be and <sup>26</sup>Al profiles show an exponential distribution with an apparent attenuation length of c. 0.5 m, similar to the spallation attenuation length in granites (ρ=2.6 g/cm<sup>2</sup>), which is Λ = 160/2.6 ≈ 60 cm. This is compatible with a near-surface accumulation, where spallation is the main process that produce cosmogenic isotopes.
- The <sup>21</sup>Ne profile also show a similar "near-surface" attenuation, excluding the top sample, that show a <sup>21</sup>Ne concentration lower than the highest sample from the core. 
- The <sup>21</sup>Ne profile is more scattered than the <sup>10</sup>Be and <sup>26</sup>Al profiles and do not converge to a 0 concentration at depth. This is a typical distribution of stable cosmonuclides profiles, as the quartz accumulates <sup>21</sup>Ne since its formation. This apparent "extra" <sup>21</sup>Ne concentration of c. 20 Matoms/g correspond, either to nucleogenic <sup>21</sup>Ne, either to <sup>21</sup>Ne accumulated several meters below the surface by muonic radiation.
- [<sup>26</sup>Al]/[<sup>10</sup>Be] ratios are between 6.3 and 10.4. This values are around surface production rate ratios (~7). These non-depleted ratios indicate that most of <sup>10</sup>Be and <sup>26</sup>Al atoms formed recently compared to their half-lifes (1.4 and 0.7 Ma respectively).
- [<sup>21</sup>Ne]/[<sup>10</sup>Be] ratios are between one and two orders of manitude higher than surface prodution rates (~4), even considering a non-cosmogenic [<sup>21</sup>Ne] of 22 Matoms/g (the lowest <sup>21</sup>Ne concentration). As the [<sup>26</sup>Al]/[<sup>10</sup>Be] ratios are not depleted and the [<sup>10</sup>Be]/[<sup>21</sup>Ne] ratios are highly depleted, we should assume that a high proportion of the <sup>21</sup>Ne atoms formed much before that most <sup>10</sup>Be and <sup>26</sup>Al atoms compared to their half-lifes.

This dataset suggests that most of the <sup>21</sup>Ne formed near the surface (near-surface attenuation distribution) before most of the <sup>10</sup>Be and <sup>26</sup>Al, and then the quartz was shielded from comic radiation enough time for most of the <sup>10</sup>Be and <sup>26</sup>Al to decay. The recent near-surface exposure formed most of the <sup>10</sup>Be and <sup>26</sup>Al (natural ratios) and a small fraction of the total <sup>21</sup>Ne concentrations.

### Surface production rates

We calculated the apparent <sup>10</sup>Be, <sup>26</sup>Al, and <sup>21</sup>Ne ages of the top sample using [The online calculators formerly known as the CRONUS-Earth online calculators v.3](https://hess.ess.washington.edu/math/v3/v3_age_in.html). We used the LSDn ages and the spreadsheet from the [Average cosmogenic production rate calculator](https://github.com/angelrodes/average_cosmogenic_production_rate_calculator) to calculte the surface production rates for <sup>10</sup>Be, <sup>26</sup>Al, and <sup>21</sup>Ne by spallation and muons, and their correspondig attenuation lengths, according to the approximations shown in [Rodés (2021)](https://doi.org/10.3390/geosciences11090362).

## Digital elevation model

As the samples are taken from inside a odd-shaped landform and the samples are exposed to cosmic rays that cross the tor in high variaty of angles and rock thicknesses, a digital elevation model (DEM) is needed to calculate how much cosmic radiation the samples received during the geological history of the tor.

### Core location

The GPS data taken in the field during the sampling campaign locates the core at the gographical coordinate 57.109150 N, -3.462608 E. However, according to [Google Maps](https://maps.app.goo.gl/z3deD5SJwoSBNCmq9) and the digital terrain model (DMT) used here, the top of the tor is located 9.7 m to the north of this location. Therefore, we used the coordinates [ 311510 , 802980.7 ] in the British National Grid for the location of the samples.

![coordinates_correction](https://github.com/angelrodes/Ben-Avon/assets/53089531/d318bd7f-49ba-442c-b67b-8e7072e1aeb3)\
_Original GPS coordinates on Google Earth and distance to the top of the tor._

![coordinates_correction_DEM](https://github.com/angelrodes/Ben-Avon/assets/53089531/b7bad8fb-0f06-447a-a1b6-5c781e9e20db)\
_Original and corrected coordinates on the DMT used._

### Sources of elevation data

In order to simulate the cosmic irradiaton of the core samples, we need a DEM of the terrain at the area where the tor is located, but also a very detailed surface model of the tor. To achieve this, we combined the following sources of terrain data:

- 5-m Digital Model of the Terrain from [Digimap](https://digimap.edina.ac.uk/).
- A 6 meter wide E-W cross section measured in the field during the sampling campaing.
- A picture showing the S-N contour of the tor. This contour was digitalised using g3data software ([Frantz, 2000](https://manpages.ubuntu.com/manpages/jammy/man1/g3data.1.html)).

![DSC_0136_g3data_S-N](https://github.com/angelrodes/Ben-Avon/assets/53089531/0fa6dd2c-26c4-4b48-a261-534f4f6220df)\
_Digitalised picture using G3data Graph Analyser. The person at the right of the picture was used as scale._

![Elevation_data_sources](https://github.com/angelrodes/Ben-Avon/assets/53089531/e10d690e-736e-477f-9f4a-bd79f41bcc8b)\
_DMT from [Digimap](https://digimap.edina.ac.uk/) and digitalised data from field notes and one picture. Sample positions are depicted by red stars._

### Radial DEM

Using all the topographic data, a polar mesh centered at the core position was generated. 100 azimuths are equally spaced every 3.6 degreees. Radial distances are logarithmic-like spaced to get a good resolution of the samples' shielding.

![3-D_DEM](https://github.com/angelrodes/Ben-Avon/assets/53089531/d0edf8c9-098f-4983-add8-4e2303495450)\
_Radial DEM generated form the three sources of topographic data. The left image shows the full DEM and the right image is restricted to the tor area. Sample locations, partially hidded by the grid, are depicted by red stars._

### Tor areas

Tor areas were defined by areas protruding more than 0.5 above a smmothed surface intepolated from areas far from the center of the grid.

![DEM_and_satellite](https://github.com/angelrodes/Ben-Avon/assets/53089531/55693274-2504-478d-9250-f2f3a9c529b7)\
_Defined tor areas compared to the satellite image from the same area._

## Erosion models

Three erosion models were build in order to test three different scenarios for the formation of the tor:

- Model A: This model considers the tor as an antecedent landform. The relief migrates down with time with a constat erosion rate for all points of the DEM. In this model, the tor is located in an area that was originally more prominent that the sourrounding areas. The surface position is calculated by adding vertical lowerings to the current DEM.
- Model B: This model simulates that the the tor is "exhumed" by differential erosion from an original smooth surface. In this model, the tor formed due to its differential erosionability (e.g. harder granite or less joints). The surface position is calculated by adding vertical lowerings to ta smooth version of the current DEM excluding the tors. In the defined tor areas, the surface is calculated as the maximum between the smooothed DEM and the current DEM.
- Model C: This model assumes that the tor is formed by lateral erosion. In this model, the tor is conserved due to its distance to the areas where ledge migration initiated. The surface positions are calculated by selecting the elevations of the current DEM and of an elevated smoothed DEM covering the tor. The selection is based in the distance to the center of the tor and the "horizontal lowering" for each evolution stage.

![ModelA](https://github.com/angelrodes/Ben-Avon/assets/53089531/279f01f0-ed9c-4cc6-9875-0961158865cc)\
_Model A: Antecedent landform._

![ModelB](https://github.com/angelrodes/Ben-Avon/assets/53089531/626c3d7b-61ec-4bc9-a453-c92bf708d0a0)\
_Model B: Exhumed landform._

![ModelC](https://github.com/angelrodes/Ben-Avon/assets/53089531/e477cf28-28f9-4cff-8306-32b47755a3bf)\
_Model C: Ledge migration._

## Shielding factors

## Glacial model

# Results

## Erosion Models

## Glacial models

# Discussion

# Conclusions

- Fitted erosion models suggest that the tor fromed by differential erosion rather than by epigenesis (inherited relief) or lateral erosion. This erosion model could correspond to:
  - Differential erosion due to compositional inhomogeneities in the granite.
  - Differential erosion due to rheological inhomogeneities (e.g. joint spacing).
  - Erosion of ice or sediments covering the tor surface (exhumation).
- Fitted glacial model indicate that the tor was exposed near the surface for at least 2 Ma and then covered by ice during part of the Pliocene and most of the Quaternary. However, the tor is in an area that has been ice-free at least since the Yonger Dryas.
- Two scenarios could explain the conditions needed to achieve long periods of shielding since the Pliocene in a place that is now ice-free and exposed to cosmic radiation:
  - The tor was at the bed of a glacial valley at some point during the Pliocene and was covered by a thick till deposit that diasappeared during the Quaternary due to lateral migration of the glacial valley.
  - The tor formed at a higher elevation during the Pliocene and covered by ice through most Plio-Quaternary glaciations while the whole area sank due to glacial isostatic adjustment.
