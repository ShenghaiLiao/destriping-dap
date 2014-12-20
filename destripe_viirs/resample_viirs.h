

void resample_viirs(float ** img, float ** lat, float ** lon, int nx, int ny, float tmin, float tmax)
{
  // nx = width of image (should be 3200 for VIIRS)
  // ny = height of image ( 5408 or 5392 for ~10 min VIIRS granule)
  // img[0..ny][0..nx]  - original image (brightness temperature)
  // lat[0..ny][0..nx]  - original latitude
  // lon[0..ny][0..nx]  - original longitude
  // ( tmin , tmax )    - range of "valid" data, subject ot resampling/interpolation
  // TODO: resample img, lat, lon 
  //
  return;
};
