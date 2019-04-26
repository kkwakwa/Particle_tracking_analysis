# Particle Tracking & Spot Analysis package

There are a couple of different things I am going to be trying do with this package

## Particle Tracking

1. Take ThunderSTORM files and convert them to a form that Trackmate will recognise

1. import Track data from either
  * Trackmate
  * Alan's Bayesian single particle software
  * Alan's Matlabtrack implimentation
  * Trackpy

  This data will be used to:

    1. calculate important statistics on spots, individual tracks and track aggregates
      * #### per- track
        * Mean displacement
        * Mean Square Deviation
        * Goodness of fit for MSD

      * ### per-segment
        * x-displacement
        * y-displacement
        * r-displacement

      * ### Per-spot
        * Intensity

    1. Identify merge points from Trackmate/bayesian  data

    2. Plot all of this data

    3. Bonus: plot individual track data as an animation and then overlay it against the original spot image data

## Spot Analysis
1. Import a stack of single molecule bleaching images, turn them into a maximum intensity projection, identify all the bleach locations and then save off the bleach data in raw camera ADU values and output a CSV with all the bleach data in it, identifying each track by its position
