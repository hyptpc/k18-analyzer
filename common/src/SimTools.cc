// SimTools.cc
#include "SimTools.hh"
#include "ReslMapMan.hh"
#include "DetectorList.hh"
#include "Hodo2Hit.hh"
#include "CDCWireMapMan.hh"
#include "HistTools.hh"

namespace{
  const DetectorList&  dlist = DetectorList::GetInstance();
  const CDCWireMapMan& gCDC  = CDCWireMapMan::GetInstance();
}
using namespace sim;
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void sim::Convert(DetectorData* detData, CDCAnalyzer* cdc)
{
  DCAnalyzer* dc=0;
  HodoAnalyzer *hodo=0;
  Convert(detData,cdc,dc,hodo);
}

void sim::Convert(DetectorData* detData, CDCAnalyzer* cdcAna, DCAnalyzer* dcAna, HodoAnalyzer* hodoAna)
{
   for(  int i=0; i<detData->detectorHitSize(); i++ ){
     DetectorHit *mchit = detData-> detectorHit(i);
     int cid       = mchit->detectorID();
     int pdg_id    = mchit->pdg();
     int parent_id = mchit->parentID();
     //     std::cout<<i<<"  "<<cid<<"  "<<layer<<std::endl;
     if( pdg_id==321 && parent_id==0 ) mchit->setTime(mchit->time()*-1);
     if( dlist.IsCounter(cid) && mchit->de() > 0.1){
       Hodo2Hit *hp = new Hodo2Hit( mchit );
       if( !hp ) continue;
       if( hodoAna && hp->IsCalculated())
	 hodoAna->AddHit(cid,hp);
       else
	 delete hp;
     }else if( dlist.IsChamber(cid) ){
       if( pdg_id==22 || pdg_id==2112 ) continue;
       DCHit *hit = new DCHit( mchit );
       if( !hit ) continue;
       int layer     = mchit->layerID();
       if( cid==DetIdCDC){
	 if(cdcAna ){	 
	   cdcAna->GetDCHC(layer).push_back(hit);
	 }else{
	   delete hit;
	 }
       }else{
	 if(dcAna){
	   // BLDC
	   dcAna->GetDCHC(cid,layer).push_back(hit);
	 }else{
	   delete hit;
	 }
       }
     }
   }
#if 0 
   {
     TString name="CDC";
     double dlbins[3]={2000,-5,15};
     for( int layer=0; layer<15; ++layer ){
       const double nw=gCDC.nwires(layer);    
       double mulbins[3]={nw+1,-0.5,nw+0.5};       
       const DCHitContainer &cont = cdcAna->GetDCHC(layer);
       int mul=cont.size();
       hist::H1(Form("%s_Mul_layer%d",name.Data(),layer),mul,mulbins);
       for(int ihit=0;ihit<mul;ihit++){
	 DCHit* hit=cont[ihit];      	
	 int wire =hit->GetWire();
	 double dl=hit->GetDriftLength();
	 hist::H1(Form("%s_HitPat_layer%d",name.Data(),layer),wire,mulbins);
	 hist::H1(Form("%s_dl_layer%d",name.Data(),layer),dl,dlbins);
       }
     }
   }
#endif   
}

