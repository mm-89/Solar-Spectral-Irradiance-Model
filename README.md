# Solar-Spectral-Irradiance-Model
This is a simple solar spectral irradiance model to simulate direct (soon diffuse) irradiance at Earth surface.

This model comes from Bird and Riordan, 1985, and it is going to use in parallel to InExEs model (https://github.com/mm-89/InExEs).

# How to use it

```bash
# upload the library
from spectralIrradiance import SpectralIrradiance
```
```bash
# create an object with default value of Ozone
# precipitable water and aerosol optical depth
current_irradiance = SpectralIrradiance()
```

```bash
# give spectral irradiance at given SZA (degree), standard
# pressure and day of the year
my.get_irradiance(45, 1020, 180)
```

```bash
# or plot
my.plot_irradiance(45, 1020, 180)
```
