#include <ConfMan.hh>
extern int
hddaq_main();
extern int
unidaq_main();
extern int
knucl_main();

namespace
{
  ConfMan&              gConf     = ConfMan::GetInstance();
}

//______________________________________________________________________________
int
main( int argc, char **argv )
{
  gConf.ParseCommand(argc,argv);
  if(gConf.IsMC())
    return knucl_main();
  else
    return hddaq_main();

  //  return EXIT_SUCCESS;
}
