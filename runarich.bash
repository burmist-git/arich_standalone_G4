make clean; make -j8; 
#./arich run.mac 123123 arich_diffPoll.root opticalphoton 0 0 -0.2 420 0 0 2
#./arich run.mac 123123 arich_randomPoll.root opticalphoton 0 0 -0.2 420 0 0 2
#./arich run.mac 123123 arich_randomPoll_250nm.root opticalphoton 0 0 -0.2 250 0 0 2
#./arich run.mac 123123 arich_randomPoll_300nm.root opticalphoton 0 0 -0.2 300 0 0 2
#./arich run.mac 123123 arich_randomPoll_500nm.root opticalphoton 0 0 -0.2 500 0 0 2
#./arich run.mac 123123 arich_randomPoll_600nm.root opticalphoton 0 0 -0.2 600 0 0 2
#./arich run.mac 123123 arich_randomPoll_650nm.root opticalphoton 0 0 -0.2 650 0 0 2
#./arich run.mac 123123 arich_mu_2cm.root mu- 0 0 -0.2 3 0 0 2
#./arich vis.mac 123123 arichV.root mu- 0 0 -1.1 3 0 0 2
#./arich vis.mac 123123 arichV.root opticalphoton 0 0 -0.2 420 0 0 2

#geomID=0
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

#geomID=1
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

#geomID=2
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

#geomID=3
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

#geomID=4
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

#geomID=5
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

#geomID=6
#echo "geomID = $geomID"
#./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
#mv G4Data0.heprep G4Data0_geomID_$geomID.heprep

geomID=7
echo "geomID = $geomID"
./arich vis.mac 123123 arichV.root mu- 0 0 0 3 33 82 $geomID
mv G4Data0.heprep G4Data0_geomID_$geomID.heprep
