/// \file
/// \ingroup tutorial_geom
/// Script drawing a detector geometry (here ATLAS).
///
/// by default the geometry is drawn using the GL viewer
/// Using the TBrowser, you can select other components
/// if the file containing the geometry is not found in the local
/// directory, it is automatically read from the ROOT web site.
///
/// \macro_code
///
/// \author Rene Brun

void geomID7_VGM_vis() {
  //TGeoManager::Import("http://root.cern.ch/files/atlas.root");
  TGeoManager::Import("geometry.root");
  //gGeoManager->DefaultColors();
  gGeoManager->SetMaxVisNodes(5000);
  //gGeoManager->SetVisLevel(4);
  gGeoManager->GetVolume("World")->Draw("ogl");
  //new TBrowser;
}
