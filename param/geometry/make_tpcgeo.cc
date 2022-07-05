void make_tpcgeo(){

  TGeoManager::Import("FC.gdml");
  TGeoVolume *fieldcageV = (TGeoVolume*) gGeoManager -> GetVolume("FcLV") -> Clone("fc");
  TGeoMaterial *fieldcageMat = (TGeoMaterial*) gGeoManager -> GetVolume("FcLV") -> GetMaterial() -> Clone();
  TGeoManager::Import("Target.gdml");
  TGeoVolume *targetV = (TGeoVolume*) gGeoManager -> GetVolume("TargetLV") -> Clone("target");
  TGeoMaterial *targetMat = (TGeoMaterial*) gGeoManager -> GetVolume("TargetLV") -> GetMaterial() -> Clone();
  TGeoManager::Import("Tpc.gdml");
  TGeoVolume *tpcV = (TGeoVolume*) gGeoManager -> GetVolume("TpcLV")-> Clone("tpc");
  TGeoMaterial *tpcMat = (TGeoMaterial*) gGeoManager -> GetVolume("TpcLV") -> GetMaterial() -> Clone();

  targetV -> SetLineColor(kBlue);
  fieldcageV -> SetLineColor(kRed);
  tpcV -> SetLineColor(kGreen);

  //Transition
  TGeoRotation *rot = new TGeoRotation();
  rot -> RotateX(-90.);
  TGeoTranslation *trans = new TGeoTranslation(0,0,-143); //cm

  TGeoVolume* hyptpc = new TGeoVolumeAssembly("hyptpc");
  hyptpc -> AddNode(tpcV,1,rot);
  hyptpc -> AddNode(targetV,1,trans);
  hyptpc -> AddNode(fieldcageV,1,rot);

  gGeoManager -> SetTopVolume(hyptpc);
  //gGeoManager -> AddMaterial(tpcMat);
  gGeoManager -> AddMaterial(targetMat);
  gGeoManager -> AddMaterial(fieldcageMat);
  gGeoManager -> SetTopVisible(1);
  gGeoManager -> CloseGeometry();
  //gGeoManager -> Export("temp_hyptpcGeo.gdml");
  gGeoManager -> Export("temp_hyptpcGeo.root");
  gGeoManager -> GetTopVolume() -> Draw("ogl");

}
