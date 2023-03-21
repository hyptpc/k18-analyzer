#ifndef TPC_VERTEX_HH
#define TPC_VERTEX_HH

/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <GFRaveVertex.h>

#include <TObject.h>

class TPCVertex {

public:

  TPCVertex();
  TPCVertex(const genfit::GFRaveVertex &vertex);
  virtual ~TPCVertex();

  void Clear();
  
  void SetIsCollisionVertex(bool flag) {isCollisionVertex = flag;}
  bool IsCollisionVertex() const {return isCollisionVertex;}
  
  void SetIsDecayVertex(bool flag) {isDecayVertex = flag;}
  bool IsDecayVertex() const {return isDecayVertex;}

  void SetVertexID(Int_t id) {fVertexID = id;}
  Int_t GetVertexID() const {return fVertexID;} 
  
  void SetPosition(TVector3 pos) {fPosition = pos;}
  TVector3 GetPosition() const {return fPosition;}

  void SetChisqr(Double_t chisqr) {fChisqr = chisqr;}
  Double_t GetChisqr() const {return fChisqr;}

  //  void SetCovariance(Double_t xCov, Double_t yCov, Double_t zCov);
  TMatrixDSym GetCovariance() const {return fCovariance;}
  
  void SetNumTrack(Int_t n) { fNumTrack = n;}
  Int_t GetNumTrack() const {return fNumTrack;}
  
private:

  bool isCollisionVertex;
  bool isDecayVertex;

  Int_t fVertexID;
  Int_t fNumTrack;

  TVector3 fPosition;
  TMatrixDSym fCovariance;
  Double_t fChisqr;

};

#endif
