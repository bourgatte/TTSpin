#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";

if($UserID eq "vcherepa"){
    $UserIDCern="cherepan";
#    $UserDir="--vcherepa";
}

$PWD=getcwd;

#printf("\n ---> Your user ID is:   $UserID \n");
if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\n\n\n ========================================================================================");
    printf("\nWelcome to tauola installer, please look at the instruction below.");
    printf("\nThis code requires one input option. The syntax is: ./todo.pl [OPTION]");
    printf("\n First of all setup the environment: \n\n");
    printf("\n./todo.pl --setup  <tauoladir>                              For example: ./todo.pl --setup  install\n");
    printf("\n\nAfter this step is completed prcoceed further and  ");
    printf("\nchoose from the following options:\n");
    printf("\n./todo.pl --help                                             Prints this message");
    printf("\n./todo.pl --tauola  <tauoladir>                              Install tauola and user codes; <tauoladir> must be the same as was set with setup option ");
    printf("\n./todo.pl --tauoladefault  <tauoladir>                       Install default  tauola; <tauoladir> must be the same as was set with setup option  ");
    printf("\n  ========================================================================================\n");
    exit(0);  
}
my $dir = getcwd;
$time= strftime("%h_%d_%Y",localtime);



for($l=0;$l<$numArgs; $l++){
    
    if($ARGV[$l] eq "--setup"){
	$setdir=$ARGV[l+1];


	system(sprintf("rm Install_TauolaEnvironment_*"));
	$SLDP='\$LD_LIBRARY_PATH';


	system(sprintf("rm Install_TauolaEnvironment_*"));	
	system(sprintf("echo \"export PYTHIA8DATA='$PWD/$setdir/tauola++/1.1.5/pythia8/176/xmldoc'\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=$SLDP:$PWD/$setdir/tauola++/1.1.5/TauSpiner/lib\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=$SLDP:$PWD/$setdir/tauola++/1.1.5/pythia8/176/lib/\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=$SLDP:$PWD/$setdir/tauola++/1.1.5/HepMC-2.06.05/workdir/lib \">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=$SLDP:$PWD/$setdir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"export LD_LIBRARY_PATH=$SLDP:$PWD/$setdir/tauola++/1.1.5/examples/UserCodes\">> Install_TauolaEnvironment_$time"));
#	system(sprintf("echo \"cernlib-use --version 5.34.18 root \n\">> Install_TauolaEnvironment_$time"));
	system(sprintf("echo \"source /libcern/root/5.34.18/sl6.3-x86_64/setup.sh \n\">> Install_TauolaEnvironment_$time"));




	system(sprintf("cp Makefile.template Makefile; "));
	system(sprintf("./subs '{DIR}'  $PWD/$setdir/  Makefile; "));

	printf("\n\nInstructions:");
	printf("\nTo complete this step do:  \n\n");
	printf("\n    1) source  Install_TauolaEnvironment_$time \n");
	printf("\n    2)  ./todo.pl --tauola $setdir \n\n");
    }
    if($ARGV[$l] eq "--tauola"){
	$tauoladir=$ARGV[l+1];
#	$l++;
	$currentdir=getcwd;
	$s1_par='\$(pwd)';
	$SP='\$outfile';

	system(sprintf("cernlib-use --version 5.34.18 root \n"));

	system(sprintf("rm -rf $tauoladir \n"));
	printf("\nInstalling Tauola++  to  $tauoladir \n");
	system(sprintf("mkdir  $tauoladir \n"));
	system(sprintf("cd $tauoladir; wget http://service-spi.web.cern.ch/service-spi/external/MCGenerators/distribution/tauola++/tauola++-1.1.5-src.tgz; tar -xzvf tauola++-1.1.5-src.tgz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://service-spi.web.cern.ch/service-spi/external/MCGenerators/distribution/pythia8/pythia8-176-src.tgz ; tar -xzvf pythia8-176-src.tgz ;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://www.hepforge.org/archive/lhapdf/lhapdf-5.9.1.tar.gz; tar -xzvf lhapdf-5.9.1.tar.gz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://mc-tester.web.cern.ch/MC-TESTER/MC-TESTER-1.25.0.tar.gz; tar -xzvf MC-TESTER-1.25.0.tar.gz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://lcgapp.cern.ch/project/simu/HepMC/download/HepMC-2.06.05.tar.gz;  tar -xzvf  HepMC-2.06.05.tar.gz;"));

	printf("\n Downloading Complete ... \n");
	printf("\n Start Installation ... \n\n\n");
	printf("\n ___________________Installing HEPMC ... _____________________\n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; .././configure -prefix=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir  -with-momentum=GEV -with-length=CM; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; make; make check; make install;"));
	printf("\n___________________Installing  LHAPDF ... _____________________\n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/; ./configure -prefix=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir  --libdir=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/; make;  make install;"));
	system(sprintf("mkdir  $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/share/lhapdf/PDFsets;"));
	system(sprintf("cd  $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/share/lhapdf/PDFsets; ../../../bin/lhapdf-getdata --repo=http://www.hepforge.org/archive/lhapdf/pdfsets/5.9.1 MSTW2008nnlo90cl.LHgrid; "));
	printf("\n ___________________Installing pythia8 ... _____________________ \n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/;  ./configure --enable-shared --lcgplatform=slc6_amd64_gcc530-opt --with-hepmc=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir --with-hepmcversion=2.06.05;"));	
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/;  make; "));
	printf("\n ___________________Compiling MC-TESTER ... _____________________\n\n\n");
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/;  ./configure --with-HepMC=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/;  make; "));
	printf("\n ___________________Compiling tauola ... _____________________\n\n\n");

	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export PYTHIA8DATA=$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/xmldoc;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauSpiner/tauola++/1.1.5/lib;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/lib/;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/;export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir/lib; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib"));
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/;   ./configure --prefix=$PWD/$tauoladir/tauola++/1.1.5/workdir  --with-hepmc=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir  --with-pythia8=$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/  --with-lhapdf=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/ --with-mc-tester=$PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/   --with-tau-spinner; "));


#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make;"));
#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make all;"));
#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make install;"));
#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples; make"));

	printf("\n___________________Tauola is compiled _____________________\n");
	printf("|                                                                                               |\n");
	printf("|                                                                                               |\n");
	printf("|                                                                                               |\n");
	printf("|                                                                                               |\n");
	printf("\n ___________________User Codes_______ _____________________\n\n\n");
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples;  mkdir UserCodes;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples;  mkdir output;"));


	system(sprintf("cp UserCodes/Makefile.shared $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/main.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/a1Helper.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/a1Helper.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/rhoHelper.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/rhoHelper.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/TauPolInterface.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/TauPolInterface.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/TauDecaysHelper.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/TauDecaysHelper.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/PolarimetricA1.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/PolarimetricA1.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));

	system(sprintf("cp UserCodes/SCalculator.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/SCalculator.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));

	system(sprintf("cp UserCodes/TauTauSpiner.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
	system(sprintf("cp UserCodes/TauTauSpiner.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));


	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes/; make  -f Makefile.shared;"));
#	system(sprintf("cp Makefile  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
#	system(sprintf("cp poltaumain_pythia_tauola.cxx  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
##	system(sprintf("cp poltaumain_pythia_tauola.Po  $PWD/$tauoladir/tauola++/1.1.5/examples/.deps/;"));
#	system(sprintf("echo \"\#dummy\">> $PWD/$tauoladir/tauola++/1.1.5/examples/.deps/poltaumain_pythia_tauola.Po"));
#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples; make"));






	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make all;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make install;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples; make"));

	system(sprintf("cp Makefile  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
	system(sprintf("cp transspin_pythia_tauola.cxx  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
	system(sprintf("cp transspin_pythia_tauola.Po  $PWD/$tauoladir/tauola++/1.1.5/examples/.deps/;"));





#	system(sprintf("echo \"\#dummy\">> $PWD/$tauoladir/tauola++/1.1.5/examples/.deps/poltaumain_pythia_tauola.Po"));



	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples; make"));





       # Setting up qsub submitter
        system(sprintf("cp submit $PWD/$tauoladir/tauola++/1.1.5/examples/; "));
        system(sprintf("cp Submit.sh $PWD/$tauoladir/tauola++/1.1.5/examples/; "));
	system(sprintf("echo \"#! /bin/bash     \">> qsub_submit.sh"));
	system(sprintf("echo \"echo 'Starting Job'       \">> qsub_submit.sh"));
	system(sprintf("echo \"export workdir=$s1_par     \">> qsub_submit.sh"));
	system(sprintf("echo \"export HOME=$s1_par     \">> qsub_submit.sh"));
	system(sprintf("echo \"cd $PWD/$tauoladir/tauola++/1.1.5/examples/;    \">> qsub_submit.sh"));
	system(sprintf("echo \"source $PWD/Install_TauolaEnvironment_$time   \">> qsub_submit.sh"));
	system(sprintf("echo \"$PWD/$tauoladir/tauola++/1.1.5/examples/poltaumain_pythia_tauola.exe $SP  \">> qsub_submit.sh"));
	system(sprintf("echo \"echo 'Completed Job'    \">> qsub_submit.sh"));
	system(sprintf("mv qsub_submit.sh $PWD/$tauoladir/tauola++/1.1.5/examples/; "));


	printf("\nInstruction:   \n    ");
	printf(" cd $PWD/$tauoladir/tauola++/1.1.5/examples/; ./transspin_pythia_tauola.exe <OutputFileName> default_output | tee log.out;   \n\n\n");
    }

    

    if($ARGV[$l] eq "--tauoladefault"){
	$tauoladir=$ARGV[l+1];
#	$l++;
	$currentdir=getcwd;


	system(sprintf("cernlib-use --version 5.34.18 root \n"));
                                

	system(sprintf("rm -rf $tauoladir \n"));
	printf("\nInstalling Tauola++  to  $tauoladir \n");
	system(sprintf("mkdir  $tauoladir \n"));
	system(sprintf("cd $tauoladir; wget http://service-spi.web.cern.ch/service-spi/external/MCGenerators/distribution/tauola++/tauola++-1.1.5-src.tgz; tar -xzvf tauola++-1.1.5-src.tgz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://service-spi.web.cern.ch/service-spi/external/MCGenerators/distribution/pythia8/pythia8-176-src.tgz ; tar -xzvf pythia8-176-src.tgz ;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://www.hepforge.org/archive/lhapdf/lhapdf-5.9.1.tar.gz; tar -xzvf lhapdf-5.9.1.tar.gz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://mc-tester.web.cern.ch/MC-TESTER/MC-TESTER-1.25.0.tar.gz; tar -xzvf MC-TESTER-1.25.0.tar.gz;"));
	system(sprintf("cd $tauoladir/tauola++/1.1.5/; wget http://lcgapp.cern.ch/project/simu/HepMC/download/HepMC-2.06.05.tar.gz;  tar -xzvf  HepMC-2.06.05.tar.gz;"));
	printf("\n Downloading Complete ... \n");
	printf("\n Start Installation ... \n\n\n");
	printf("\n ___________________Installing HEPMC ... _____________________\n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; .././configure -prefix=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir  -with-momentum=GEV -with-length=CM; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; make; make check; make install;"));
	printf("\n___________________Installing  LHAPDF ... _____________________\n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/; ./configure -prefix=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir  --libdir=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/; make;  make install;"));
	system(sprintf("mkdir  $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/share/lhapdf/PDFsets;"));
	system(sprintf("cd  $PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/share/lhapdf/PDFsets; ../../../bin/lhapdf-getdata --repo=http://www.hepforge.org/archive/lhapdf/pdfsets/5.9.1 MSTW2008nnlo90cl.LHgrid; "));
	printf("\n ___________________Installing pythia8 ... _____________________ \n\n\n");
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/;  ./configure --enable-shared --lcgplatform=slc6_amd64_gcc530-opt --with-hepmc=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir --with-hepmcversion=2.06.05;"));	
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/pythia8/176/;  make; "));
	printf("\n ___________________Compiling MC-TESTER ... _____________________\n\n\n");
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/;  ./configure --with-HepMC=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/;  make; "));
	printf("\n ___________________Compiling tauola ... _____________________\n\n\n");

	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export PYTHIA8DATA=$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/xmldoc;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauSpiner/tauola++/1.1.5/lib;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/lib/;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/;export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir/lib; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/lib"));
	system(sprintf("mkdir $PWD/$tauoladir/tauola++/1.1.5/workdir; "));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/;   ./configure --prefix=$PWD/$tauoladir/tauola++/1.1.5/workdir  --with-hepmc=$PWD/$tauoladir/tauola++/1.1.5/HepMC-2.06.05/workdir  --with-pythia8=$PWD/$tauoladir/tauola++/1.1.5/pythia8/176/  --with-lhapdf=$PWD/$tauoladir/tauola++/1.1.5/lhapdf-5.9.1/workdir/ --with-mc-tester=$PWD/$tauoladir/tauola++/1.1.5/MC-TESTER/   --with-tau-spinner; "));
#	system(sprintf("cp Makefile  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
#	system(sprintf("cp mypythia_example.cxx  $PWD/$tauoladir/tauola++/1.1.5/examples;"));
#	system(sprintf("echo \"\#dummy\">> $PWD/$tauoladir/tauola++/1.1.5/examples/.deps/mypythia_example.Po"));
#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples;  mkdir UserCodes;"));
#	system(sprintf("cp UserCodes/Makefile.shared $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
#	system(sprintf("cp UserCodes/MultiplyNumbers.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
#	system(sprintf("cp UserCodes/MultiplyNumbers.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
#	system(sprintf("cp UserCodes/AddNumbers.h $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
#	system(sprintf("cp UserCodes/AddNumbers.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
#	system(sprintf("cp UserCodes/main.cc $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes; "));
#	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples/UserCodes/; make all -f Makefile.shared;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make all;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/; make install;"));
	system(sprintf("cd $PWD/$tauoladir/tauola++/1.1.5/examples; make"));


	printf("\nInstruction:   \n    ");
	printf(" cd $PWD/$tauoladir/tauola++/1.1.5/examples/; ./taumain_pythia_example | tee log.out;   \n\n\n");
    }




}

