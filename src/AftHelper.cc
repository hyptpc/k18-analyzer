// -*- C++ -*-

#include "AftHelper.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <cmath>
#include <stdint.h>
#include <cstdlib>

#include "FuncName.hh"

// const Double_t TpcPadHelper::PadParameter[32][6]=
// { // Layer#, #ofPad, division#, radius#, padLength
//   {  0,  48,  14.75,  48, 0,  9.0 },
//   {  1,  48,  24.25,  48, 0,  9.0 },
//   {  2,  72,  33.75,  72, 0,  9.0 },
//   {  3,  96,  43.25,  96, 0,  9.0 },
//   {  4, 120,  52.75, 120, 0,  9.0 },
//   {  5, 144,  62.25, 144, 0,  9.0 },
//   {  6, 168,  71.75, 168, 0,  9.0 },
//   {  7, 192,  81.25, 192, 0,  9.0 },
//   {  8, 216,  90.75, 216, 0,  9.0 },
//   {  9, 240, 100.25, 240, 0,  9.0 },
//   { 10, 208, 111.50, 241, 0, 12.5 },
//   { 11, 218, 124.50, 271, 0, 12.5 },
//   { 12, 230, 137.50, 300, 0, 12.5 },
//   { 13, 214, 150.50, 330, 0, 12.5 },
//   { 14, 212, 163.50, 360, 0, 12.5 },
//   { 15, 214, 176.50, 390, 0, 12.5 },
//   { 16, 220, 189.50, 420, 0, 12.5 },
//   { 17, 224, 202.50, 449, 0, 12.5 },
//   { 18, 232, 215.50, 479, 0, 12.5 },
//   { 19, 238, 228.50, 509, 0, 12.5 },
//   { 20, 244, 241.50, 539, 0, 12.5 },
//   { 21, 232, 254.50, 569, 0, 12.5 },
//   { 22, 218, 267.50, 599, 0, 12.5 },
//   { 23, 210, 280.50, 628, 0, 12.5 },
//   { 24, 206, 293.50, 658, 0, 12.5 },
//   { 25, 202, 306.50, 688, 0, 12.5 },
//   { 26, 200, 319.50, 718, 0, 12.5 },
//   { 27, 196, 332.50, 748, 0, 12.5 },
//   { 28, 178, 345.50, 777, 0, 12.5 },
//   { 29, 130, 358.50, 807, 0, 12.5 },
//   { 30, 108, 371.50, 837, 0, 12.5 },
//   { 31,  90, 384.50, 867, 0, 12.5 } };

// //_____________________________________________________________________________
// TpcPadParam::TpcPadParam( Int_t asad, Int_t aget, Int_t channel,
//                           Int_t pad_id, Int_t layer_id, Int_t row_id )
//   : m_asad( asad ),
//     m_aget( aget ),
//     m_channel( channel ),
//     m_pad_id( pad_id ),
//     m_layer_id( layer_id ),
//     m_row_id( row_id ),
//     m_key( Key( asad, aget, channel ) )
// {
// }

// //_____________________________________________________________________________
// TpcPadParam::~TpcPadParam( void )
// {
// }

// //_____________________________________________________________________________
// void
// TpcPadParam::Print( void ) const
// {
//   std::cout << " AsAd=" << m_asad
// 	    << " AGet=" << m_aget
// 	    << " Channel=" << m_channel
// 	    << " PadId=" << m_pad_id
// 	    << " LayerId=" << m_layer_id
// 	    << " RowId=" << m_row_id
// 	    << " Key=" << m_key
// 	    << std::endl;
// }


//_____________________________________________________________________________
AftHelper::AftHelper( void )
{
}

//_____________________________________________________________________________
AftHelper::~AftHelper( void )
{
}

// //_____________________________________________________________________________
// Int_t
// TpcPadHelper::GetPadId( Int_t asad, Int_t aget, Int_t channel ) const
// {
//   auto param = GetParam( asad, aget, channel );
//   if( param && param->AsAdId()>=0 )
//     return param->PadId();
//   else
//     return -1;
// }

// //_____________________________________________________________________________
// Int_t
// TpcPadHelper::GetPadId( Int_t layer, Int_t row ) const
// {
//   auto param = GetParam( layer, row );
//   if( param && param->AsAdId()>=0 )
//     return param->PadId();
//   else
//     return -1;
// }

// //_____________________________________________________________________________
// TpcPadParam*
// TpcPadHelper::GetParam( Int_t asad, Int_t aget, Int_t channel ) const
// {
//   Int_t key = TpcPadParam::Key( asad, aget, channel );
//   TpcPadParam* param = nullptr;
//   auto itr = m_map.find( key );
//   if( itr != m_map.end() )
//     param = itr->second;
// #if 0
//   else
//     std::cerr << FUNC_NAME << " No parameter, AsAd=" << asad
// 	      << ", AGet=" << aget << ", Channel=" << channel << std::endl;
// #endif
//   return param;
// }

// //_____________________________________________________________________________
// TpcPadParam*
// TpcPadHelper::GetParam( Int_t layer, Int_t row ) const
// {
//   for( const auto& p : m_map ){
//     if( p.second->LayerId() == layer &&
//         p.second->RowId() == row ){
//       return p.second;
//     }
//   }
// #if 0
//   std::cerr << FUNC_NAME << " No parameter, Layer=" << layer
//             << ", Row=" << row << std::endl;
// #endif
//   return nullptr;
// }

// //_____________________________________________________________________________
// TpcPadParam*
// TpcPadHelper::GetParam( Int_t pad ) const
// {
//   for( const auto& p : m_map ){
//     if( p.second->PadId() == pad )
//       return p.second;
//   }
// #if 0
//   std::cerr << FUNC_NAME << " No parameter, Pad=" << pad << std::endl;
// #endif
//   return nullptr;
// }

//_____________________________________________________________________________
Double_t
AftHelper::GetX( Int_t plane, Int_t seg )
{
  Double_t centerx = m_posSeg0[plane%4] + seg*m_distance_seg;
  return centerx;
}

//_____________________________________________________________________________
Double_t
AftHelper::GetZ( Int_t plane, Int_t seg )
{

  Double_t centerz = m_posZ0;
  centerz += (plane+1)/2*m_distance_x_xp + plane/2*m_distance_xp_y;
  return centerz;
}

// //_____________________________________________________________________________
// TVector3
// TpcPadHelper::GetPoint( Int_t pad ) const
// {
//   TVector3 vec;
//   Int_t layer, row;
//   Int_t sum=0;
//   for( layer=0; layer<=30 && sum+TpcPadHelper::PadParameter[layer][1]<=pad; layer++ ){
//     sum += TpcPadHelper::PadParameter[layer][1];
//   }
//   row = pad - sum;

//   if( row>TpcPadHelper::PadParameter[layer][1] ) {
//     vec.SetXYZ(0,-1,0);
//   } else {
//     Double_t stheta = 180.-(360./TpcPadHelper::PadParameter[layer][3])*TpcPadHelper::PadParameter[layer][1]/2.;
//     Double_t theta = stheta+(row+0.5)*360./TpcPadHelper::PadParameter[layer][3]-180.;

//     Double_t x = TpcPadHelper::PadParameter[layer][2]*sin(theta*acos(-1)/180.);
//     Double_t z = TpcPadHelper::PadParameter[layer][2]*cos(theta*acos(-1)/180.)-143.;
//     vec.SetXYZ(x,0,z);
//   }
//   return vec;
// }

// //_____________________________________________________________________________
// Bool_t
// TpcPadHelper::Initialize( const TString &file_name )
// {
//   std::ifstream ifs( file_name );
//   if( !ifs.is_open() ){
//     std::cerr << FUNC_NAME << " " << file_name
// 	      << " ... failed to open" << std::endl;
//     return false;
//   }
//   TString line;
//   while( ifs.good() && line.ReadLine( ifs ) ){
//     if( line[0]=='#' )
//       continue;
//     std::istringstream iss( line.Data() );
//     Int_t aget, asad, channel, layer, row, pad;
//     iss >> aget >> asad >> channel >> layer >> row >> pad;
//     if( aget < 0 || asad < 0 || channel < 0 )
//       continue;
//     Int_t key = TpcPadParam::Key( asad, aget, channel );
//     // std::cout << line << std::endl;
//     if( m_map.find( key ) != m_map.end() ){
//       std::cerr << FUNC_NAME << " " << file_name
// 		<< " ... duplicated parameter" << std::endl;
//       TpcPadParam( asad, aget, channel, pad, layer, row ).Print();
//       m_map.find( key )->second->Print();
//       continue;
//       // return false;
//     }
//     if( GetParam( pad ) ){
//       std::cerr << FUNC_NAME << " " << file_name
// 		<< " ... duplicated parameter" << std::endl;
//       TpcPadParam( asad, aget, channel, pad, layer, row ).Print();
//       GetParam( pad )->Print();
//       continue;
//       // return false;
//     }
//     m_map[key] = new TpcPadParam( asad, aget, channel,
//                                   pad, layer, row );
//   }
// #if 0
//   std::cout << FUNC_NAME << " " << file_name
// 	    << " ... initialized (NParam=" << m_map.size() << ")" << std::endl;
//   for( const auto& p : m_map ){
//     std::cout << " PadId=" << p.second->PadId()
// 	      << " LayerId=" << p.second->LayerId()
// 	      << " RowId=" << p.second->RowId()
// 	      << std::endl;
//   }
// #endif
//   return true;
// }
