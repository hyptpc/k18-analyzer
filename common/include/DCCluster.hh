/* DCCluster.hh */
#ifndef DCCluster_h
#define DCCluster_h 1
#include <iostream>
#include <vector>
#include <TVector3.h>

class DCHit;

class DCCluster
{ 
private:
  typedef std::vector<DCHit* > DCHitContainer;
  DCHitContainer m_hit;
  std::vector<int> m_nthhit;

public:
  DCCluster();
  virtual ~DCCluster();
  
 private:
  //  int m_cluster_id;
  double m_time;
  double m_timesub;
  double m_ctime;
  TVector3 m_pos;
  TVector3 m_dir;

 public:
  void Calc( const bool &isMC=false );
 
  // int GetClusterID() const { return m_cluster_id; }
  // int id() const { return m_cluster_id; }
  double GetTime() const { return m_time; }
  double time() const { return m_time; }
  double GetTimeSub() const { return m_timesub; }
  double timesub() const { return m_timesub; }
  double GetCTime() const { return m_ctime; }
  double ctime() const { return m_ctime; }
  double x() const { return m_pos.X(); } 
  double y() const { return m_pos.Y(); } 
  double z() const { return m_pos.Z(); } 
  TVector3 pos() const { return m_pos; }
  TVector3 dir() const { return m_dir; }

  //  void SetClusterID( const short &i ){ ClusterID = i; }
  DCHit* hit(int nh) { return m_hit[nh]; }
  int nth(int nh) { return m_nthhit[nh]; }
  void SetHit( DCHit* hit, int nth ){ m_hit.push_back(hit); m_nthhit.push_back(nth); }
  int  nhit() const { return m_hit.size(); }
  void Clear(); 
};
#endif
