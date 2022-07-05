void monitor_tpcgeo(){

  TGeoManager::Import("hyptpcGeo.gdml");
  gGeoManager -> GetTopVolume() -> Draw("ogl");

}
